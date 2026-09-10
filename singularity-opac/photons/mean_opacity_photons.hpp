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
#ifndef SINGULARITY_OPAC_PHOTONS_MEAN_OPACITY_PHOTONS_
#define SINGULARITY_OPAC_PHOTONS_MEAN_OPACITY_PHOTONS_
// This file was made in part with generative AI.

#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>
#include <singularity-opac/base/opac_error.hpp>
#include <singularity-opac/base/robust_utils.hpp>
#include <singularity-opac/base/sp5.hpp>
#include <singularity-opac/constants/constants.hpp>
#include <singularity-opac/photons/mean_photon_types.hpp>
#include <singularity-opac/photons/mean_photon_utils.hpp>
#include <singularity-opac/photons/mean_photon_variant.hpp>
#include <singularity-opac/photons/non_cgs_photons.hpp>
#include <singularity-opac/photons/thermal_distributions_photons.hpp>
#include <spiner/databox.hpp>

namespace singularity {
namespace photons {
namespace impl {

// TODO(BRR) Note: It is assumed that lambda is constant for all densities and
// temperatures

template <typename pc = PhysicalConstantsCGS>
class MeanOpacity {
 public:
  using PC = pc;
  using DataBox = Spiner::DataBox<Real>;

  MeanOpacity() = default;

  template <typename Opacity, typename GroupBoundsIndexer>
  MeanOpacity(const Opacity &opac, const Real lRhoMin, const Real lRhoMax,
              const int NRho, const Real lTMin, const Real lTMax, const int NT,
              const GroupBoundsIndexer &group_bounds, const int ngroups,
              const int NNuPerGroup = 64, Real *lambda = nullptr) {
    MeanOpacityImpl_(opac, lRhoMin, lRhoMax, NRho, lTMin, lTMax, NT,
                     group_bounds, ngroups, NNuPerGroup, lambda);
  }

  template <typename GroupBoundsIndexer>
  MeanOpacity(const DataBox &kappaPlanck, const DataBox &kappaRosseland,
              const GroupBoundsIndexer &group_bounds) {
    // Table-backed multigroup opacities always carry explicit group bounds.
    // To represent [nu_max, infinity), the final bound must literally be
    // IEEE +infinity, not a large finite proxy value.
    LoadOpacityTables_(kappaPlanck.size() > 0 ? &kappaPlanck : nullptr,
                       kappaRosseland.size() > 0 ? &kappaRosseland : nullptr,
                       group_bounds);
  }

  template <typename GroupBoundsIndexer>
  MeanOpacity(const DataBox &opacity, const int gmode,
              const GroupBoundsIndexer &group_bounds) {
    if (gmode != Planck && gmode != Rosseland) {
      OPAC_ERROR("photons::MeanOpacity: invalid opacity averaging mode");
    }
    LoadOpacityTables_(gmode == Planck ? &opacity : nullptr,
                       gmode == Rosseland ? &opacity : nullptr, group_bounds);
  }

#ifdef SPINER_USE_HDF
  MeanOpacity(const std::string &filename, const int matid) {
    LoadHDF_(filename, matid);
  }

  MeanOpacity(const std::string &filename, const std::string &material_name) {
    LoadHDF_(filename, material_name);
  }

  void Save(const std::string &filename, const std::string &material_name,
            const bool append = false) const {
    if (material_name.empty()) {
      OPAC_ERROR("photons::MeanOpacity: material name must not be empty");
    }
    Save_(filename, "/" + material_name, nullptr, material_name, append);
  }

  void Save(const std::string &filename, const int matid,
            const bool append = false) const {
    Save(filename, matid, std::string(), append);
  }

  void Save(const std::string &filename, const int matid,
            const std::string &material_name,
            const bool append = false) const {
    Save_(filename, "/" + std::to_string(matid), &matid, material_name,
          append);
  }

  void Save(const std::string &filename, const int matid,
            const char *material_name, const bool append = false) const {
    Save(filename, matid, std::string(material_name), append);
  }

 private:
  void Save_(const std::string &filename, const std::string &material_path,
             const int *matid, const std::string &material_name,
             const bool append) const {
    DataBox kappaPlanck;
    DataBox kappaRosseland;
    DataBox groupBounds;
    ExportOpacityTables_(kappaPlanck, kappaRosseland);
    ExportGroupBounds(groupBounds, groupBounds_, ngroups_);

    herr_t status = H5_SUCCESS;
    hid_t file = append ? H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT)
                        : H5Fcreate(filename.c_str(), H5F_ACC_TRUNC,
                                    H5P_DEFAULT, H5P_DEFAULT);
    hid_t material = -1;
    if (append && H5Lexists(file, material_path.c_str(), H5P_DEFAULT) > 0) {
      material = H5Gopen(file, material_path.c_str(), H5P_DEFAULT);
    } else {
      material = H5Gcreate(file, material_path.c_str(), H5P_DEFAULT,
                           H5P_DEFAULT, H5P_DEFAULT);
    }
    if (matid != nullptr) {
      status += H5LTset_attribute_int(file, material_path.c_str(),
                                      SP5::Material::matid, matid, 1);
    }
    if (!material_name.empty()) {
      status += H5LTset_attribute_string(file, material_path.c_str(),
                                          SP5::Material::name,
                                          material_name.c_str());
    }
    if (hasPlanck_) {
      status += kappaPlanck.saveHDF(
          material, SP5::MultigroupOpac::PlanckGroupOpacity);
    }
    if (hasRosseland_) {
      status += kappaRosseland.saveHDF(
          material, SP5::MultigroupOpac::RosselandGroupOpacity);
    }
    // Absorption and scattering share group bounds when appended.
    if (H5Lexists(material, SP5::Multigroup::GroupBounds, H5P_DEFAULT) >
        0) {
      DataBox existingBounds;
      const herr_t bounds_status =
          existingBounds.loadHDF(material, SP5::Multigroup::GroupBounds);
      if (bounds_status != H5_SUCCESS ||
          existingBounds.size() != groupBounds.size()) {
        existingBounds.finalize();
        OPAC_ERROR("photons::MeanOpacity: existing material group bounds "
                   "are incompatible with appended opacity tables");
      }
      for (int i = 0; i < groupBounds.size(); ++i) {
        if (existingBounds(i) != groupBounds(i)) {
          existingBounds.finalize();
          OPAC_ERROR("photons::MeanOpacity: existing material group bounds "
                     "are incompatible with appended opacity tables");
        }
      }
      existingBounds.finalize();
    } else {
      status += groupBounds.saveHDF(material, SP5::Multigroup::GroupBounds);
    }
    status += H5Gclose(material);
    status += H5Fclose(file);

    if (hasPlanck_) kappaPlanck.finalize();
    if (hasRosseland_) kappaRosseland.finalize();
    groupBounds.finalize();

    if (status != H5_SUCCESS) {
      OPAC_ERROR("photons::MeanOpacity: HDF5 error\n");
    }
  }

 public:
#endif

  PORTABLE_INLINE_FUNCTION
  void PrintParams() const {
    printf("Photon multigroup opacity. ngroups = %d\n", ngroups_);
  }

  MeanOpacity GetOnDevice() {
    MeanOpacity other;
    if (hasPlanck_) {
      other.lkappaPlanck_ = Spiner::getOnDeviceDataBox(lkappaPlanck_);
    }
    if (hasRosseland_) {
      other.lkappaRosseland_ = Spiner::getOnDeviceDataBox(lkappaRosseland_);
    }
    other.groupBounds_ = Spiner::getOnDeviceDataBox(groupBounds_);
    other.ngroups_ = ngroups_;
    other.hasPlanck_ = hasPlanck_;
    other.hasRosseland_ = hasRosseland_;
    return other;
  }

  void Finalize() {
    if (hasPlanck_) lkappaPlanck_.finalize();
    if (hasRosseland_) lkappaRosseland_.finalize();
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

  PORTABLE_INLINE_FUNCTION
  bool HasPlanckOpacity() const noexcept { return hasPlanck_; }

  PORTABLE_INLINE_FUNCTION
  bool HasRosselandOpacity() const noexcept { return hasRosseland_; }

  std::vector<Real> GetGroupBounds() const {
    std::vector<Real> bounds(ngroups_ + 1);
    for (int group = 0; group <= ngroups_; ++group) {
      bounds[group] = groupBounds_(group);
    }
    return bounds;
  }

  // The group-index-less "mean" accessors only make sense when there is a
  // single group. In that case the lone group spans the entire spectrum, so
  // its group-integrated coefficient (stored at group 0) IS the traditional
  // gray mean. We therefore require ngroups==1 and forward to group 0. This is
  // not a distinguished "mean slot": for ngroups>1, group 0 is simply the
  // lowest-frequency group and callers must use the group-index API.
  PORTABLE_INLINE_FUNCTION
  Real PlanckMeanAbsorptionCoefficient(const Real rho, const Real temp) const {
    PORTABLE_REQUIRE(
        ngroups_ == 1,
        "PlanckMeanAbsorptionCoefficient only valid for ngroups==1. "
        "Use PlanckGroupAbsorptionCoefficient(rho, temp, group) for "
        "multigroup.");
    return PlanckGroupAbsorptionCoefficient(rho, temp, 0);
  }

  // See PlanckMeanAbsorptionCoefficient: the gray mean is the single group 0.
  PORTABLE_INLINE_FUNCTION
  Real RosselandMeanAbsorptionCoefficient(const Real rho,
                                          const Real temp) const {
    PORTABLE_REQUIRE(
        ngroups_ == 1,
        "RosselandMeanAbsorptionCoefficient only valid for ngroups==1. "
        "Use RosselandGroupAbsorptionCoefficient(rho, temp, group) for "
        "multigroup.");
    return RosselandGroupAbsorptionCoefficient(rho, temp, 0);
  }

  PORTABLE_INLINE_FUNCTION
  Real PlanckGroupAbsorptionCoefficient(const Real rho, const Real temp,
                                        const int group) const {
    PORTABLE_REQUIRE(hasPlanck_,
                     "photons::MeanOpacity: Planck opacity is unavailable");
    return GroupAbsorptionCoefficient_(lkappaPlanck_, rho, temp, group);
  }

  PORTABLE_INLINE_FUNCTION
  Real RosselandGroupAbsorptionCoefficient(const Real rho, const Real temp,
                                           const int group) const {
    PORTABLE_REQUIRE(
        hasRosseland_,
        "photons::MeanOpacity: Rosseland opacity is unavailable");
    return GroupAbsorptionCoefficient_(lkappaRosseland_, rho, temp, group);
  }

  PORTABLE_INLINE_FUNCTION
  Real AbsorptionCoefficient(const Real rho, const Real temp, const int group,
                             const int gmode = Rosseland) const {
    return (gmode == Planck)
               ? PlanckGroupAbsorptionCoefficient(rho, temp, group)
               : RosselandGroupAbsorptionCoefficient(rho, temp, group);
  }

  // Logarithmic temperature derivative of the group absorption coefficient,
  // d(log alpha_g)/d(log T) at fixed rho, where alpha_g = rho * kappa_g
  PORTABLE_INLINE_FUNCTION
  Real PlanckGroupDLogAbsorptionCoefficientDLogT(const Real rho,
                                                 const Real temp,
                                                 const int group) const {
    PORTABLE_REQUIRE(hasPlanck_,
                     "photons::MeanOpacity: Planck opacity is unavailable");
    return GroupDLogAbsCoeffDLogT_(lkappaPlanck_, rho, temp, group);
  }

  PORTABLE_INLINE_FUNCTION
  Real RosselandGroupDLogAbsorptionCoefficientDLogT(const Real rho,
                                                    const Real temp,
                                                    const int group) const {
    PORTABLE_REQUIRE(
        hasRosseland_,
        "photons::MeanOpacity: Rosseland opacity is unavailable");
    return GroupDLogAbsCoeffDLogT_(lkappaRosseland_, rho, temp, group);
  }

  PORTABLE_INLINE_FUNCTION
  Real DLogAbsorptionCoefficientDLogT(const Real rho, const Real temp,
                                      const int group,
                                      const int gmode = Rosseland) const {
    return (gmode == Planck)
               ? PlanckGroupDLogAbsorptionCoefficientDLogT(rho, temp, group)
               : RosselandGroupDLogAbsorptionCoefficientDLogT(rho, temp, group);
  }

  // Like the mean accessors above, Emissivity has no group-index argument and
  // so is only defined for ngroups==1, where group 0 is the whole-spectrum
  // (gray) group.
  PORTABLE_INLINE_FUNCTION
  Real Emissivity(const Real rho, const Real temp, const int gmode = Rosseland,
                  Real *lambda = nullptr) const {
    if (ngroups_ != 1) {
      OPAC_ERROR("photons::MeanOpacity: Emissivity only valid for ngroups==1");
    }
    PlanckDistribution<PC> dist;
    Real B = dist.ThermalDistributionOfT(temp, lambda);
    return AbsorptionCoefficient(rho, temp, 0, gmode) * B;
  }

  PORTABLE_INLINE_FUNCTION
  int GroupOfNu(const Real nu) const {
    if (!(nu >= GroupBoundAt(groupBounds_, 0) &&
          nu <= GroupBoundAt(groupBounds_, ngroups_))) {
      OPAC_ERROR("photons::MeanOpacity: frequency is outside group bounds");
    }
    return GroupOfNuImpl(groupBounds_, ngroups_, nu);
  }

  PORTABLE_INLINE_FUNCTION
  Real PlanckGroupAbsorptionCoefficientFromNu(const Real rho, const Real temp,
                                              const Real nu) const {
    return AbsorptionCoefficientFromNu(rho, temp, nu, Planck);
  }

  PORTABLE_INLINE_FUNCTION
  Real RosselandGroupAbsorptionCoefficientFromNu(const Real rho,
                                                 const Real temp,
                                                 const Real nu) const {
    return AbsorptionCoefficientFromNu(rho, temp, nu, Rosseland);
  }

  PORTABLE_INLINE_FUNCTION
  Real AbsorptionCoefficientFromNu(const Real rho, const Real temp,
                                   const Real nu,
                                   const int gmode = Rosseland) const {
    return AbsorptionCoefficient(rho, temp, GroupOfNu(nu), gmode);
  }

 private:
#ifdef SPINER_USE_HDF
  void LoadHDF_(const std::string &filename, const int matid) {
    LoadHDF_(filename, OpenMaterialGroupByMatid, matid);
  }

  void LoadHDF_(const std::string &filename,
                const std::string &material_name) {
    LoadHDF_(filename, OpenMaterialGroupByName, material_name);
  }

  template <typename MaterialOpener, typename Selector>
  void LoadHDF_(const std::string &filename, MaterialOpener opener,
                const Selector &selector) {
    DataBox kappaPlanck;
    DataBox kappaRosseland;
    DataBox groupBounds;
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
      OPAC_ERROR("photons::MeanOpacity: unable to open HDF5 file");
    }
    hid_t material = opener(file, selector);
    if (material < 0) {
      OPAC_ERROR("photons::MeanOpacity: material group not found in HDF5 file");
    }

    bool has_planck = LoadOpacityDataBoxIfPresent_(
        material, SP5::MultigroupOpac::PlanckGroupOpacity, kappaPlanck);
    bool has_rosseland = LoadOpacityDataBoxIfPresent_(
        material, SP5::MultigroupOpac::RosselandGroupOpacity, kappaRosseland);
    bool has_gray_rosseland = false;
    if (!has_rosseland) {
      // A gray Rosseland absorption table is a valid one-group table.  It is
      // promoted to rank 3 so the runtime has one uniform representation.
      has_rosseland = LoadOpacityDataBoxIfPresent_(
          material, SP5::MeanOpac::RosselandMeanOpacity, kappaRosseland);
      has_gray_rosseland = has_rosseland;
    }
    const herr_t bounds_status =
        LoadGroupBounds(material, SP5::Multigroup::GroupBounds,
                        groupBounds);
    H5Gclose(material);
    H5Fclose(file);

    if (bounds_status != H5_SUCCESS || (!has_planck && !has_rosseland)) {
      OPAC_ERROR("photons::MeanOpacity: HDF5 error\n");
    }

    if (has_gray_rosseland && groupBounds.size() != 2) {
      groupBounds.finalize();
      groupBounds.resize(2);
      groupBounds(0) = 0.;
      groupBounds(1) = std::numeric_limits<Real>::infinity();
    }

    LoadOpacityTables_(has_planck ? &kappaPlanck : nullptr,
                       has_rosseland ? &kappaRosseland : nullptr,
                       groupBounds);
    groupBounds.finalize();
    if (has_planck) kappaPlanck.finalize();
    if (has_rosseland) kappaRosseland.finalize();
  }

  bool LoadOpacityDataBoxIfPresent_(const hid_t material, const char *field,
                                    DataBox &data) const {
    if (H5Lexists(material, field, H5P_DEFAULT) <= 0) return false;
    DataBox source;
    if (source.loadHDF(material, field) != H5_SUCCESS) return false;
    if (source.rank() == 3) {
      data.copyMetadata(source);
      for (int i = 0; i < source.size(); ++i) data(i) = source(i);
    } else if (source.rank() == 2) {
      data.resize(source.dim(2), source.dim(1), 1);
      data.setRange(1, source.range(0));
      data.setRange(2, source.range(1));
      for (int iRho = 0; iRho < source.dim(2); ++iRho) {
        for (int iT = 0; iT < source.dim(1); ++iT) {
          data(iRho, iT, 0) = source(iRho, iT);
        }
      }
    } else {
      source.finalize();
      return false;
    }
    source.finalize();
    return true;
  }
#endif

  PORTABLE_INLINE_FUNCTION
  Real GroupAbsorptionCoefficient_(const DataBox &lkappa, const Real rho,
                                   const Real temp, const int group) const {
    const Real lRho = ToLog(rho);
    const Real lT = ToLog(temp);
    return rho * FromLog(lkappa.interpToReal(lRho, lT, group));
  }

  // Evaluate slope d(log kappa)/d(log T) via table lookups, which equals
  // d(log alpha)/d(log T) (the log-base and the rho prefactor drop out of the
  // ratio).
  PORTABLE_INLINE_FUNCTION
  Real GroupDLogAbsCoeffDLogT_(const DataBox &lkappa, const Real rho,
                               const Real temp, const int group) const {
    const Real lRho = ToLog(rho);
    const Real lT = ToLog(temp);
    const auto lT_grid =
        lkappa.range(1); // axis 1 is log10(T) (see setRange in Impl_)
    const Real dlT = lT_grid.dx();
    const int iT = lT_grid.index(lT); // clamped to a valid cell [iT, iT+1]
    const Real lT_lo = lT_grid.x(iT);
    const Real lT_hi = lT_grid.x(iT + 1);
    const Real L_lo = lkappa.interpToReal(lRho, lT_lo, group);
    const Real L_hi = lkappa.interpToReal(lRho, lT_hi, group);
    return (L_hi - L_lo) / dlT;
  }

  void ValidateOpacityTable_(const DataBox &table,
                             const char *averaging_name) const {
    if (table.rank() != 3) {
      OPAC_ERROR("photons::MeanOpacity: opacity tables must be rank 3");
    }
    if (table.dim(1) <= 0) {
      OPAC_ERROR("photons::MeanOpacity: ngroups must be positive");
    }
    (void)averaging_name;
  }

  void ValidateOpacityCompatibility_(const DataBox &first,
                                     const DataBox &second) const {
    for (int dim = 1; dim <= 3; ++dim) {
      if (first.dim(dim) != second.dim(dim)) {
        OPAC_ERROR("photons::MeanOpacity: table dimensions do not match");
      }
    }
    if (first.range(1) != second.range(1) ||
        first.range(2) != second.range(2)) {
      OPAC_ERROR("photons::MeanOpacity: table ranges do not match");
    }
  }

  template <typename GroupBoundsIndexer>
  void LoadOpacityTables_(const DataBox *kappaPlanck,
                          const DataBox *kappaRosseland,
                          const GroupBoundsIndexer &group_bounds) {
    if (kappaPlanck == nullptr && kappaRosseland == nullptr) {
      OPAC_ERROR("photons::MeanOpacity: at least one opacity table is required");
    }
    if (kappaPlanck != nullptr) ValidateOpacityTable_(*kappaPlanck, "Planck");
    if (kappaRosseland != nullptr)
      ValidateOpacityTable_(*kappaRosseland, "Rosseland");
    if (kappaPlanck != nullptr && kappaRosseland != nullptr) {
      ValidateOpacityCompatibility_(*kappaPlanck, *kappaRosseland);
    }
    const DataBox &reference =
        kappaPlanck != nullptr ? *kappaPlanck : *kappaRosseland;
    ngroups_ = reference.dim(1);
    ValidateGroupBounds(group_bounds, ngroups_);
    SetGroupBounds(groupBounds_, group_bounds, ngroups_);
    hasPlanck_ = kappaPlanck != nullptr;
    hasRosseland_ = kappaRosseland != nullptr;
    if (hasPlanck_) {
      lkappaPlanck_.copyMetadata(*kappaPlanck);
      for (int i = 0; i < kappaPlanck->size(); ++i) {
        lkappaPlanck_(i) = ToLog((*kappaPlanck)(i));
      }
    }
    if (hasRosseland_) {
      lkappaRosseland_.copyMetadata(*kappaRosseland);
      for (int i = 0; i < kappaRosseland->size(); ++i) {
        lkappaRosseland_(i) = ToLog((*kappaRosseland)(i));
      }
    }
  }

  void ExportOpacityTables_(DataBox &kappaPlanck,
                            DataBox &kappaRosseland) const {
    if (hasPlanck_) {
      kappaPlanck.copyMetadata(lkappaPlanck_);
      for (int i = 0; i < lkappaPlanck_.size(); ++i) {
        kappaPlanck(i) = FromLog(lkappaPlanck_(i));
      }
    }
    if (hasRosseland_) {
      kappaRosseland.copyMetadata(lkappaRosseland_);
      for (int i = 0; i < lkappaRosseland_.size(); ++i) {
        kappaRosseland(i) = FromLog(lkappaRosseland_(i));
      }
    }
  }

  template <typename Opacity, typename GroupBoundsIndexer>
  void MeanOpacityImpl_(const Opacity &opac, const Real lRhoMin,
                        const Real lRhoMax, const int NRho, const Real lTMin,
                        const Real lTMax, const int NT,
                        const GroupBoundsIndexer &group_bounds,
                        const int ngroups, const int NNuPerGroup,
                        Real *lambda = nullptr) {
#ifndef NDEBUG
    auto RPC = RuntimePhysicalConstants(PC());
    auto opc = opac.GetRuntimePhysicalConstants();
    assert(RPC == opc && "Physical constants are the same");
#endif

    if (NNuPerGroup < 2) {
      OPAC_ERROR("photons::MeanOpacity: NNuPerGroup must be at least 2");
    }
    ValidateGroupBounds(group_bounds, ngroups);

    ngroups_ = ngroups;
    hasPlanck_ = true;
    hasRosseland_ = true;
    SetGroupBounds(groupBounds_, group_bounds, ngroups_);
    lkappaPlanck_.resize(NRho, NT, ngroups_);
    lkappaPlanck_.setRange(1, lTMin, lTMax, NT);
    lkappaPlanck_.setRange(2, lRhoMin, lRhoMax, NRho);
    lkappaRosseland_.copyMetadata(lkappaPlanck_);

    PlanckDistribution<PC> dist;
    std::vector<Real> planckDenom(ngroups_, 0.);
    std::vector<Real> rosselandDenom(ngroups_, 0.);

    for (int iT = 0; iT < NT; ++iT) {
      const Real lT = lkappaPlanck_.range(1).x(iT);
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
        const Real lRho = lkappaPlanck_.range(2).x(iRho);
        const Real rho = FromLog(lRho);

        for (int group = 0; group < ngroups_; ++group) {
          Real kappaPlanckNum = 0.;
          Real kappaRosselandNum = 0.;
          const Real nuMin = GroupBoundAt(group_bounds, group);
          const Real nuMax = GroupBoundAt(group_bounds, group + 1);
          ForEachGroupFrequencySample<PC>(
              T, nuMin, nuMax, NNuPerGroup, [&](const Real nu, const Real dnu) {
                const Real alpha =
                    opac.AbsorptionCoefficient(rho, T, nu, lambda);
                Real B = 0.;
                Real dBdT = 0.;
                ThermalWeightsAtNu<PC>(dist, T, nu, B, dBdT);
                kappaPlanckNum += alpha / rho * B * dnu;

                if (alpha > singularity_opac::robust::SMALL()) {
                  kappaRosselandNum +=
                      singularity_opac::robust::ratio(rho, alpha) * dBdT * dnu;
                }
              });

          const Real kappaPlanck = singularity_opac::robust::ratio(
              kappaPlanckNum, planckDenom[group]);
          const Real kappaRosseland =
              (rosselandDenom[group] > singularity_opac::robust::SMALL() &&
               kappaRosselandNum > singularity_opac::robust::SMALL())
                  ? singularity_opac::robust::ratio(rosselandDenom[group],
                                                    kappaRosselandNum)
                  : 0.;

          lkappaPlanck_(iRho, iT, group) = ToLog(kappaPlanck);
          lkappaRosseland_(iRho, iT, group) = ToLog(kappaRosseland);
          if (std::isnan(lkappaPlanck_(iRho, iT, group)) ||
              std::isnan(lkappaRosseland_(iRho, iT, group))) {
            OPAC_ERROR("photons::MeanOpacity: NAN in opacity evaluations");
          }
        }
      }
    }
  }

  DataBox lkappaPlanck_;
  DataBox lkappaRosseland_;
  DataBox groupBounds_;
  int ngroups_ = 0;
  bool hasPlanck_ = false;
  bool hasRosseland_ = false;
};

} // namespace impl

using MeanOpacityBase = impl::MeanOpacity<PhysicalConstantsCGS>;
using MeanOpacity =
    impl::MeanVariant<MeanOpacityBase, MeanNonCGSUnits<MeanOpacityBase>>;

} // namespace photons
} // namespace singularity

#endif // SINGULARITY_OPAC_PHOTONS_MEAN_OPACITY_PHOTONS_
