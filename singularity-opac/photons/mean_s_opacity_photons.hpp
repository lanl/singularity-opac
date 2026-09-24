// ======================================================================
// © 2022-2026. Triad National Security, LLC. All rights reserved.  This
// program was produced under U.S. Government contract
// 89233218CNA000001 for Los Alamos National Laboratory (LANL), which
// is operated by Triad National Security, LLC for the U.S.
// Department of Energy/National Nuclear Security Administration. All
// rights in the program are reserved by Triad National Security, LLC,
// and the U.S. Department of Energy/National Nuclear Security
// Administration. The Government is granted for itself and others
// acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
// license in this material to reproduce, prepare derivative works,
// distribute copies to the public, perform publicly and display
// publicly, and to permit others to do so.
// ======================================================================
#ifndef SINGULARITY_OPAC_PHOTONS_MEAN_S_OPACITY_PHOTONS_
#define SINGULARITY_OPAC_PHOTONS_MEAN_S_OPACITY_PHOTONS_
// This file was made in part with generative AI.

#include <cassert>
#include <cmath>
#include <cstdio>
#include <optional>
#include <string>
#include <vector>

#include <ports-of-call/portability.hpp>
#include <singularity-opac/base/opac_error.hpp>
#include <singularity-opac/base/robust_utils.hpp>
#include <singularity-opac/base/sp5.hpp>
#include <singularity-opac/constants/constants.hpp>
#include <singularity-opac/photons/mean_photon_types.hpp>
#include <singularity-opac/photons/mean_photon_utils.hpp>
#include <singularity-opac/photons/non_cgs_s_photons.hpp>
#include <singularity-opac/photons/thermal_distributions_photons.hpp>
#include <spiner/databox.hpp>

namespace singularity {
namespace photons {
namespace impl {

// TODO(BRR) Note: It is assumed that lambda is constant for all densities and
// temperatures

template <typename pc = PhysicalConstantsCGS>
class MeanSOpacity {
 public:
  using PC = pc;
  using DataBox = Spiner::DataBox<Real>;

  MeanSOpacity() = default;

  template <typename SOpacity, typename GroupBoundsIndexer>
  MeanSOpacity(const SOpacity &s_opac, const Real lRhoMin, const Real lRhoMax,
               const int NRho, const Real lTMin, const Real lTMax, const int NT,
               const GroupBoundsIndexer &group_bounds, const int ngroups,
               const int NNuPerGroup = 64, Real *lambda = nullptr) {
    MeanSOpacityImpl_(s_opac, lRhoMin, lRhoMax, NRho, lTMin, lTMax, NT,
                      group_bounds, ngroups, NNuPerGroup, lambda);
  }

  template <typename GroupBoundsIndexer>
  MeanSOpacity(const DataBox &sigmaPlanck, const DataBox &sigmaRosseland,
               const GroupBoundsIndexer &group_bounds) {
    LoadScatteringTables_(sigmaPlanck.size() > 0 ? &sigmaPlanck : nullptr,
                          sigmaRosseland.size() > 0 ? &sigmaRosseland : nullptr,
                          group_bounds);
  }

  template <typename GroupBoundsIndexer>
  MeanSOpacity(const DataBox &sigma, const int gmode,
               const GroupBoundsIndexer &group_bounds) {
    if (gmode != Planck && gmode != Rosseland) {
      OPAC_ERROR("photons::MeanSOpacity: invalid opacity averaging mode");
    }
    LoadScatteringTables_(gmode == Planck ? &sigma : nullptr,
                          gmode == Rosseland ? &sigma : nullptr, group_bounds);
  }

#ifdef SPINER_USE_HDF
  MeanSOpacity(const std::string &filename, const int opacid) {
    LoadHDF_(filename, opacid);
  }

  MeanSOpacity(const std::string &filename, const std::string &material_name) {
    LoadHDF_(filename, material_name);
  }

  void Save(const std::string &filename, const std::string &material_name,
            const bool append = false) const {
    if (material_name.empty()) {
      OPAC_ERROR("photons::MeanSOpacity: material name must not be empty");
    }
    Save_(filename, "/" + material_name, std::nullopt, material_name, append);
  }

  void Save(const std::string &filename, const int opacid,
            const bool append = false) const {
    Save(filename, opacid, std::string(), append);
  }

  // A name, when supplied, keys the material group; the opacid is recorded as
  // metadata so that several groups may share one opacid.
  void Save(const std::string &filename, const int opacid,
            const std::string &material_name, const bool append = false) const {
    const std::string material_path = material_name.empty()
                                          ? "/" + std::to_string(opacid)
                                          : "/" + material_name;
    Save_(filename, material_path, opacid, material_name, append);
  }

  void Save(const std::string &filename, const int opacid,
            const char *material_name, const bool append = false) const {
    Save(filename, opacid, std::string(material_name), append);
  }

 private:
  void Save_(const std::string &filename, const std::string &material_path,
             const std::optional<int> opacid, const std::string &material_name,
             const bool append) const {
    DataBox sigmaPlanck;
    DataBox sigmaRosseland;
    DataBox groupBounds;
    ExportScatteringTables_(sigmaPlanck, sigmaRosseland);
    ExportGroupBounds(groupBounds, groupBounds_, ngroups_);

    ScopedH5ErrorHandler h5_errors;
    hid_t file = OpenFileWrite(filename, append);
    hid_t material = CreateOrOpenGroup(file, material_path, append);
    SetAttribute(file, material_path, SP5::Offsets::opac_messageName,
                 SP5::Offsets::opac_message);
    if (opacid.has_value()) {
      SetAttribute(file, material_path, SP5::Material::opacid, *opacid);
    }
    if (!material_name.empty()) {
      SetAttribute(file, material_path, SP5::Material::opac_name,
                   material_name);
    }
    if (HasPlanckSOpacity()) {
      SaveDataBox(material, SP5::MultigroupSOpac::PlanckGroupSOpacity,
                  sigmaPlanck);
    }
    if (HasRosselandSOpacity()) {
      SaveDataBox(material, SP5::MultigroupSOpac::RosselandGroupSOpacity,
                  sigmaRosseland);
    }
    // Absorption and scattering share group bounds when appended.
    bool bounds_compatible = true;
    if (H5Lexists(material, SP5::Multigroup::GroupBounds, H5P_DEFAULT) > 0) {
      DataBox existingBounds;
      LoadDataBox(material, SP5::Multigroup::GroupBounds, existingBounds);
      bounds_compatible = existingBounds.size() == groupBounds.size();
      for (int i = 0; bounds_compatible && i < groupBounds.size(); ++i) {
        bounds_compatible = existingBounds(i) == groupBounds(i);
      }
      existingBounds.finalize();
    } else {
      SaveDataBox(material, SP5::Multigroup::GroupBounds, groupBounds);
    }
    CloseGroup(material);
    CloseFile(file);

    sigmaPlanck.finalize();
    sigmaRosseland.finalize();
    groupBounds.finalize();

    if (!bounds_compatible) {
      OPAC_ERROR("photons::MeanSOpacity: existing material group bounds are "
                 "incompatible with appended scattering tables");
    }
  }

 public:
#endif

  PORTABLE_INLINE_FUNCTION
  void PrintParams() const {
    printf("Photon multigroup scattering opacity. ngroups = %d\n", ngroups_);
  }

  MeanSOpacity GetOnDevice() {
    MeanSOpacity other;
    other.lsigmaPlanck_ = Spiner::getOnDeviceDataBox(lsigmaPlanck_);
    other.lsigmaRosseland_ = Spiner::getOnDeviceDataBox(lsigmaRosseland_);
    other.groupBounds_ = Spiner::getOnDeviceDataBox(groupBounds_);
    other.ngroups_ = ngroups_;
    return other;
  }

  void Finalize() {
    lsigmaPlanck_.finalize();
    lsigmaRosseland_.finalize();
    groupBounds_.finalize();
  }

  PORTABLE_INLINE_FUNCTION RuntimePhysicalConstants
  GetRuntimePhysicalConstants() const {
    return RuntimePhysicalConstants(PC());
  }

  PORTABLE_INLINE_FUNCTION
  int ngroups() const noexcept { return ngroups_; }

  PORTABLE_INLINE_FUNCTION
  bool HasGroupBounds() const noexcept { return true; }

  // A table that was never loaded is default constructed, and so has rank 0.
  PORTABLE_INLINE_FUNCTION
  bool HasPlanckSOpacity() const noexcept { return lsigmaPlanck_.rank() > 0; }

  PORTABLE_INLINE_FUNCTION
  bool HasRosselandSOpacity() const noexcept {
    return lsigmaRosseland_.rank() > 0;
  }

  std::vector<Real> GetGroupBounds() const {
    std::vector<Real> bounds(ngroups_ + 1);
    for (int group = 0; group <= ngroups_; ++group) {
      bounds[group] = groupBounds_(group);
    }
    return bounds;
  }

  // With ngroups==1 the lone group spans the spectrum, so group 0 is the gray
  // mean. There is no distinguished "mean slot" for ngroups>1.
  PORTABLE_INLINE_FUNCTION
  Real PlanckMeanScatteringCoefficient(const Real rho, const Real temp) const {
    PORTABLE_REQUIRE(
        ngroups_ == 1,
        "PlanckMeanScatteringCoefficient only valid for ngroups==1. "
        "Use PlanckGroupScatteringCoefficient(rho, temp, group) for "
        "multigroup.");
    return PlanckGroupScatteringCoefficient(rho, temp, 0);
  }

  // See PlanckMeanScatteringCoefficient: the gray mean is the single group 0.
  PORTABLE_INLINE_FUNCTION
  Real RosselandMeanScatteringCoefficient(const Real rho,
                                          const Real temp) const {
    PORTABLE_REQUIRE(
        ngroups_ == 1,
        "RosselandMeanScatteringCoefficient only valid for ngroups==1. "
        "Use RosselandGroupScatteringCoefficient(rho, temp, group) for "
        "multigroup.");
    return RosselandGroupScatteringCoefficient(rho, temp, 0);
  }

  PORTABLE_INLINE_FUNCTION
  Real PlanckGroupScatteringCoefficient(const Real rho, const Real temp,
                                        const int group) const {
    PORTABLE_REQUIRE(
        HasPlanckSOpacity(),
        "photons::MeanSOpacity: Planck scattering opacity is unavailable");
    return GroupScatteringCoefficient_(lsigmaPlanck_, rho, temp, group);
  }

  PORTABLE_INLINE_FUNCTION
  Real RosselandGroupScatteringCoefficient(const Real rho, const Real temp,
                                           const int group) const {
    PORTABLE_REQUIRE(
        HasRosselandSOpacity(),
        "photons::MeanSOpacity: Rosseland scattering opacity is unavailable");
    return GroupScatteringCoefficient_(lsigmaRosseland_, rho, temp, group);
  }

  PORTABLE_INLINE_FUNCTION
  Real ScatteringCoefficient(const Real rho, const Real temp, const int group,
                             const int gmode = Rosseland) const {
    return (gmode == Planck)
               ? PlanckGroupScatteringCoefficient(rho, temp, group)
               : RosselandGroupScatteringCoefficient(rho, temp, group);
  }

  // Logarithmic temperature derivative of the group scattering coefficient,
  // d(log alpha_g)/d(log T) at fixed rho, where alpha_g = rho * sigma_g
  PORTABLE_INLINE_FUNCTION
  Real PlanckGroupDLogScatteringCoefficientDLogT(const Real rho,
                                                 const Real temp,
                                                 const int group) const {
    PORTABLE_REQUIRE(
        HasPlanckSOpacity(),
        "photons::MeanSOpacity: Planck scattering opacity is unavailable");
    return GroupDLogSCoeffDLogT_(lsigmaPlanck_, rho, temp, group);
  }

  PORTABLE_INLINE_FUNCTION
  Real RosselandGroupDLogScatteringCoefficientDLogT(const Real rho,
                                                    const Real temp,
                                                    const int group) const {
    PORTABLE_REQUIRE(
        HasRosselandSOpacity(),
        "photons::MeanSOpacity: Rosseland scattering opacity is unavailable");
    return GroupDLogSCoeffDLogT_(lsigmaRosseland_, rho, temp, group);
  }

  PORTABLE_INLINE_FUNCTION
  Real DLogScatteringCoefficientDLogT(const Real rho, const Real temp,
                                      const int group,
                                      const int gmode = Rosseland) const {
    return (gmode == Planck)
               ? PlanckGroupDLogScatteringCoefficientDLogT(rho, temp, group)
               : RosselandGroupDLogScatteringCoefficientDLogT(rho, temp, group);
  }

  PORTABLE_INLINE_FUNCTION
  int GroupOfNu(const Real nu) const {
    if (!(nu >= GroupBoundAt(groupBounds_, 0) &&
          nu <= GroupBoundAt(groupBounds_, ngroups_))) {
      OPAC_ERROR("photons::MeanSOpacity: frequency is outside group bounds");
    }
    return GroupOfNuImpl(groupBounds_, ngroups_, nu);
  }

  PORTABLE_INLINE_FUNCTION
  Real PlanckGroupScatteringCoefficientFromNu(const Real rho, const Real temp,
                                              const Real nu) const {
    return ScatteringCoefficientFromNu(rho, temp, nu, Planck);
  }

  PORTABLE_INLINE_FUNCTION
  Real RosselandGroupScatteringCoefficientFromNu(const Real rho,
                                                 const Real temp,
                                                 const Real nu) const {
    return ScatteringCoefficientFromNu(rho, temp, nu, Rosseland);
  }

  PORTABLE_INLINE_FUNCTION
  Real ScatteringCoefficientFromNu(const Real rho, const Real temp,
                                   const Real nu,
                                   const int gmode = Rosseland) const {
    return ScatteringCoefficient(rho, temp, GroupOfNu(nu), gmode);
  }

 private:
#ifdef SPINER_USE_HDF
  void LoadHDF_(const std::string &filename, const int opacid) {
    MaterialSelector selector;
    selector.opacid = opacid;
    LoadHDF_(filename, selector);
  }

  void LoadHDF_(const std::string &filename, const std::string &material_name) {
    MaterialSelector selector;
    selector.name = material_name;
    LoadHDF_(filename, selector);
  }

  void LoadHDF_(const std::string &filename, const MaterialSelector &selector) {
    DataBox sigmaPlanck;
    DataBox sigmaRosseland;
    DataBox groupBounds;
    ScopedH5ErrorHandler h5_errors;
    hid_t file = OpenFileRead(filename);
    hid_t material = OpenMaterialGroup(file, selector);
    if (material < 0) {
      const bool ambiguous = material == MaterialAmbiguous;
      CloseFile(file);
      if (ambiguous) {
        OPAC_ERROR("photons::MeanSOpacity: several material groups share the "
                   "requested opacid; select by name instead");
      }
      OPAC_ERROR(
          "photons::MeanSOpacity: material group not found in HDF5 file");
    }

    const bool has_planck = LoadSOpacityDataBoxIfPresent_(
        material, SP5::MultigroupSOpac::PlanckGroupSOpacity, sigmaPlanck);
    const bool has_rosseland = LoadSOpacityDataBoxIfPresent_(
        material, SP5::MultigroupSOpac::RosselandGroupSOpacity, sigmaRosseland);
    LoadDataBox(material, SP5::Multigroup::GroupBounds, groupBounds);
    CloseGroup(material);
    CloseFile(file);

    if (!has_planck && !has_rosseland) {
      OPAC_ERROR("photons::MeanSOpacity: no scattering table found in material "
                 "group\n");
    }
    if (has_planck) ValidateScatteringTable_(sigmaPlanck);
    if (has_rosseland) ValidateScatteringTable_(sigmaRosseland);

    // ValidateGroupBounds and SetGroupBounds both walk group_bounds(0 ...
    // ngroups).
    const int file_ngroups =
        has_planck ? sigmaPlanck.dim(1) : sigmaRosseland.dim(1);
    if (groupBounds.size() != file_ngroups + 1) {
      OPAC_ERROR("photons::MeanSOpacity: group bounds count is inconsistent "
                 "with the scattering table group count");
    }

    LoadScatteringTables_(has_planck ? &sigmaPlanck : nullptr,
                          has_rosseland ? &sigmaRosseland : nullptr,
                          groupBounds);
    groupBounds.finalize();
    sigmaPlanck.finalize();
    sigmaRosseland.finalize();
  }

  bool LoadSOpacityDataBoxIfPresent_(const hid_t material, const char *field,
                                     DataBox &data) const {
    if (H5Lexists(material, field, H5P_DEFAULT) <= 0) return false;
    LoadDataBox(material, field, data);
    return true;
  }
#endif

  PORTABLE_INLINE_FUNCTION
  Real GroupScatteringCoefficient_(const DataBox &lsigma, const Real rho,
                                   const Real temp, const int group) const {
    const Real lRho = ToLog(rho);
    const Real lT = ToLog(temp);
    return rho * FromLog(lsigma.interpToReal(lRho, lT, group));
  }

  // Evaluate slope d(log sigma)/d(log T) via table lookups, which equals
  // d(log alpha)/d(log T) (the log-base and the rho prefactor drop out of the
  // ratio).
  PORTABLE_INLINE_FUNCTION
  Real GroupDLogSCoeffDLogT_(const DataBox &lsigma, const Real rho,
                             const Real temp, const int group) const {
    const Real lRho = ToLog(rho);
    const Real lT = ToLog(temp);
    const auto lT_grid =
        lsigma.range(1); // axis 1 is log10(T) (see setRange in Impl_)
    const Real dlT = lT_grid.dx();
    const int iT = lT_grid.index(lT); // clamped to a valid cell [iT, iT+1]
    const Real lT_lo = lT_grid.x(iT);
    const Real lT_hi = lT_grid.x(iT + 1);
    const Real L_lo = lsigma.interpToReal(lRho, lT_lo, group);
    const Real L_hi = lsigma.interpToReal(lRho, lT_hi, group);
    return (L_hi - L_lo) / dlT;
  }

  void ValidateScatteringTable_(const DataBox &table) const {
    if (table.rank() != 3) {
      OPAC_ERROR("photons::MeanSOpacity: scattering tables must be rank 3");
    }
    if (table.dim(1) <= 0) {
      OPAC_ERROR("photons::MeanSOpacity: ngroups must be positive");
    }
    if (table.dim(2) < 2) {
      OPAC_ERROR("photons::MeanSOpacity: scattering tables need at least two "
                 "temperature points");
    }
    if (table.dim(3) < 2) {
      OPAC_ERROR("photons::MeanSOpacity: scattering tables need at least two "
                 "density points");
    }
  }

  void ValidateScatteringCompatibility_(const DataBox &first,
                                        const DataBox &second) const {
    for (int dim = 1; dim <= 3; ++dim) {
      if (first.dim(dim) != second.dim(dim)) {
        OPAC_ERROR("photons::MeanSOpacity: table dimensions do not match");
      }
    }
    if (first.range(1) != second.range(1) ||
        first.range(2) != second.range(2)) {
      OPAC_ERROR("photons::MeanSOpacity: table ranges do not match");
    }
  }

  template <typename GroupBoundsIndexer>
  void LoadScatteringTables_(const DataBox *sigmaPlanck,
                             const DataBox *sigmaRosseland,
                             const GroupBoundsIndexer &group_bounds) {
    if (sigmaPlanck == nullptr && sigmaRosseland == nullptr) {
      OPAC_ERROR(
          "photons::MeanSOpacity: at least one scattering table is required");
    }
    if (sigmaPlanck != nullptr) ValidateScatteringTable_(*sigmaPlanck);
    if (sigmaRosseland != nullptr) ValidateScatteringTable_(*sigmaRosseland);
    if (sigmaPlanck != nullptr && sigmaRosseland != nullptr) {
      ValidateScatteringCompatibility_(*sigmaPlanck, *sigmaRosseland);
    }
    const DataBox &reference =
        sigmaPlanck != nullptr ? *sigmaPlanck : *sigmaRosseland;
    ngroups_ = reference.dim(1);
    ValidateGroupBounds(group_bounds, ngroups_);
    SetGroupBounds(groupBounds_, group_bounds, ngroups_);
    if (sigmaPlanck != nullptr) {
      lsigmaPlanck_.copyMetadata(*sigmaPlanck);
      for (int i = 0; i < sigmaPlanck->size(); ++i) {
        lsigmaPlanck_(i) = ToLog((*sigmaPlanck)(i));
      }
    }
    if (sigmaRosseland != nullptr) {
      lsigmaRosseland_.copyMetadata(*sigmaRosseland);
      for (int i = 0; i < sigmaRosseland->size(); ++i) {
        lsigmaRosseland_(i) = ToLog((*sigmaRosseland)(i));
      }
    }
  }

  void ExportScatteringTables_(DataBox &sigmaPlanck,
                               DataBox &sigmaRosseland) const {
    if (HasPlanckSOpacity()) {
      sigmaPlanck.copyMetadata(lsigmaPlanck_);
      for (int i = 0; i < lsigmaPlanck_.size(); ++i) {
        sigmaPlanck(i) = FromLog(lsigmaPlanck_(i));
      }
    }
    if (HasRosselandSOpacity()) {
      sigmaRosseland.copyMetadata(lsigmaRosseland_);
      for (int i = 0; i < lsigmaRosseland_.size(); ++i) {
        sigmaRosseland(i) = FromLog(lsigmaRosseland_(i));
      }
    }
  }

  template <typename SOpacity, typename GroupBoundsIndexer>
  void MeanSOpacityImpl_(const SOpacity &s_opac, const Real lRhoMin,
                         const Real lRhoMax, const int NRho, const Real lTMin,
                         const Real lTMax, const int NT,
                         const GroupBoundsIndexer &group_bounds,
                         const int ngroups, const int NNuPerGroup,
                         Real *lambda = nullptr) {
    if (NNuPerGroup < 2) {
      OPAC_ERROR("photons::MeanSOpacity: NNuPerGroup must be at least 2");
    }
    ValidateGroupBounds(group_bounds, ngroups);

    ngroups_ = ngroups;
    SetGroupBounds(groupBounds_, group_bounds, ngroups_);
    lsigmaPlanck_.resize(NRho, NT, ngroups_);
    lsigmaPlanck_.setRange(1, lTMin, lTMax, NT);
    lsigmaPlanck_.setRange(2, lRhoMin, lRhoMax, NRho);
    lsigmaRosseland_.copyMetadata(lsigmaPlanck_);

    PlanckDistribution<PC> dist;
    std::vector<Real> planckDenom(ngroups_, 0.);
    std::vector<Real> rosselandDenom(ngroups_, 0.);

    for (int iT = 0; iT < NT; ++iT) {
      const Real lT = lsigmaPlanck_.range(1).x(iT);
      const Real T = FromLog(lT);

      for (int group = 0; group < ngroups_; ++group) {
        Real Baccum = 0.;
        Real dBdTaccum = 0.;
        const Real nuMin = GroupBoundAt(group_bounds, group);
        const Real nuMax = GroupBoundAt(group_bounds, group + 1);
        ForEachGroupFrequencySample<PC>(
            T, nuMin, nuMax, NNuPerGroup, [&](const Real nu, const Real dnu) {
              Real B = 0.;
              Real dBdT = 0.;
              ThermalWeightsAtNu<PC>(dist, T, nu, B, dBdT);
              Baccum += B * dnu;
              dBdTaccum += dBdT * dnu;
            });

        planckDenom[group] = Baccum;
        rosselandDenom[group] = dBdTaccum;
      }

      for (int iRho = 0; iRho < NRho; ++iRho) {
        const Real lRho = lsigmaPlanck_.range(2).x(iRho);
        const Real rho = FromLog(lRho);

        for (int group = 0; group < ngroups_; ++group) {
          Real sigmaPlanckNum = 0.;
          Real sigmaRosselandNum = 0.;
          const Real nuMin = GroupBoundAt(group_bounds, group);
          const Real nuMax = GroupBoundAt(group_bounds, group + 1);
          ForEachGroupFrequencySample<PC>(
              T, nuMin, nuMax, NNuPerGroup, [&](const Real nu, const Real dnu) {
                const Real sigma =
                    s_opac.TotalScatteringCoefficient(rho, T, nu, lambda);
                Real B = 0.;
                Real dBdT = 0.;
                ThermalWeightsAtNu<PC>(dist, T, nu, B, dBdT);
                sigmaPlanckNum += sigma / rho * B * dnu;

                if (sigma > singularity_opac::robust::SMALL()) {
                  sigmaRosselandNum +=
                      singularity_opac::robust::ratio(rho, sigma) * dBdT * dnu;
                }
              });

          const Real sigmaPlanck = singularity_opac::robust::ratio(
              sigmaPlanckNum, planckDenom[group]);
          const Real sigmaRosseland =
              (rosselandDenom[group] > singularity_opac::robust::SMALL() &&
               sigmaRosselandNum > singularity_opac::robust::SMALL())
                  ? singularity_opac::robust::ratio(rosselandDenom[group],
                                                    sigmaRosselandNum)
                  : 0.;

          lsigmaPlanck_(iRho, iT, group) = ToLog(sigmaPlanck);
          lsigmaRosseland_(iRho, iT, group) = ToLog(sigmaRosseland);
          if (std::isnan(lsigmaPlanck_(iRho, iT, group)) ||
              std::isnan(lsigmaRosseland_(iRho, iT, group))) {
            OPAC_ERROR("photons::MeanSOpacity: NAN in opacity evaluations");
          }
        }
      }
    }
  }

  DataBox lsigmaPlanck_;
  DataBox lsigmaRosseland_;
  DataBox groupBounds_;
  int ngroups_ = 0;
};

} // namespace impl

using MeanSOpacityBase = impl::MeanSOpacity<PhysicalConstantsCGS>;

} // namespace photons
} // namespace singularity

#endif // SINGULARITY_OPAC_PHOTONS_MEAN_S_OPACITY_PHOTONS_
