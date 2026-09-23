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
#include <string>
#include <type_traits>
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

  template <typename OpacityDataBox, typename GroupBoundsIndexer,
            typename std::enable_if<
                std::is_same<typename std::decay<OpacityDataBox>::type,
                             DataBox>::value,
                int>::type = 0>
  MeanSOpacity(OpacityDataBox &&sigmaRosseland,
               const GroupBoundsIndexer &group_bounds) {
    LoadScatteringTables_(sigmaRosseland, group_bounds);
  }

#ifdef SPINER_USE_HDF
  MeanSOpacity(const std::string &filename, const int matid) {
    LoadHDF_(filename, matid);
  }

  MeanSOpacity(const std::string &filename, const std::string &material_name) {
    LoadHDF_(filename, material_name);
  }

  void Save(const std::string &filename, const std::string &material_name,
            const bool append = false) const {
    if (material_name.empty()) {
      OPAC_ERROR("photons::MeanSOpacity: material name must not be empty");
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
    DataBox sigmaRosseland;
    DataBox groupBounds;
    ExportScatteringTables_(sigmaRosseland);
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
                                      SP5::Material::opacid, matid, 1);
    }
    if (!material_name.empty()) {
      status += H5LTset_attribute_string(file, material_path.c_str(),
                                          SP5::Material::opac_name,
                                          material_name.c_str());
    }
    status += sigmaRosseland.saveHDF(
        material, SP5::MultigroupSOpac::RosselandGroupSOpacity);
    // Absorption and scattering share group bounds when appended.
    if (H5Lexists(material, SP5::Multigroup::GroupBounds, H5P_DEFAULT) >
        0) {
      DataBox existingBounds;
      const herr_t bounds_status =
          existingBounds.loadHDF(material, SP5::Multigroup::GroupBounds);
      if (bounds_status != H5_SUCCESS ||
          existingBounds.size() != groupBounds.size()) {
        existingBounds.finalize();
        OPAC_ERROR("photons::MeanSOpacity: existing material group bounds "
                   "are incompatible with appended scattering tables");
      }
      for (int i = 0; i < groupBounds.size(); ++i) {
        if (existingBounds(i) != groupBounds(i)) {
          existingBounds.finalize();
          OPAC_ERROR("photons::MeanSOpacity: existing material group bounds "
                     "are incompatible with appended scattering tables");
        }
      }
      existingBounds.finalize();
    } else {
      status +=
          groupBounds.saveHDF(material, SP5::Multigroup::GroupBounds);
    }
    status += H5Gclose(material);
    status += H5Fclose(file);

    sigmaRosseland.finalize();
    groupBounds.finalize();

    if (status != H5_SUCCESS) {
      OPAC_ERROR("photons::MeanSOpacity: HDF5 error\n");
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
    other.lsigmaRosseland_ = Spiner::getOnDeviceDataBox(lsigmaRosseland_);
    other.groupBounds_ = Spiner::getOnDeviceDataBox(groupBounds_);
    other.ngroups_ = ngroups_;
    return other;
  }

  void Finalize() {
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
  Real RosselandGroupScatteringCoefficient(const Real rho, const Real temp,
                                           const int group) const {
    return GroupScatteringCoefficient_(lsigmaRosseland_, rho, temp, group);
  }

  PORTABLE_INLINE_FUNCTION
  Real ScatteringCoefficient(const Real rho, const Real temp,
                             const int group) const {
    return RosselandGroupScatteringCoefficient(rho, temp, group);
  }

  // Logarithmic temperature derivative of the group scattering coefficient,
  // d(log alpha_g)/d(log T) at fixed rho, where alpha_g = rho * sigma_g
  PORTABLE_INLINE_FUNCTION
  Real RosselandGroupDLogScatteringCoefficientDLogT(const Real rho,
                                                    const Real temp,
                                                    const int group) const {
    return GroupDLogSCoeffDLogT_(lsigmaRosseland_, rho, temp, group);
  }

  PORTABLE_INLINE_FUNCTION
  Real DLogScatteringCoefficientDLogT(const Real rho, const Real temp,
                                      const int group) const {
    return RosselandGroupDLogScatteringCoefficientDLogT(rho, temp, group);
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
  Real RosselandGroupScatteringCoefficientFromNu(const Real rho,
                                                 const Real temp,
                                                 const Real nu) const {
    return ScatteringCoefficientFromNu(rho, temp, nu);
  }

  PORTABLE_INLINE_FUNCTION
  Real ScatteringCoefficientFromNu(const Real rho, const Real temp,
                                   const Real nu) const {
    return ScatteringCoefficient(rho, temp, GroupOfNu(nu));
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
    DataBox sigmaRosseland;
    DataBox groupBounds;
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
      OPAC_ERROR("photons::MeanSOpacity: unable to open HDF5 file");
    }
    hid_t material = opener(file, selector);
    if (material < 0) {
      OPAC_ERROR(
          "photons::MeanSOpacity: material group not found in HDF5 file");
    }
    herr_t status = sigmaRosseland.loadHDF(
        material, SP5::MultigroupSOpac::RosselandGroupSOpacity);
    const herr_t bounds_status =
        LoadGroupBounds(material, SP5::Multigroup::GroupBounds,
                        groupBounds);
    H5Gclose(material);
    H5Fclose(file);

    if (bounds_status != H5_SUCCESS || status != H5_SUCCESS) {
      OPAC_ERROR("photons::MeanSOpacity: HDF5 error\n");
    }

    LoadScatteringTables_(sigmaRosseland, groupBounds);
    groupBounds.finalize();
    sigmaRosseland.finalize();
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

  void ValidateScatteringTable_(const DataBox &sigmaRosseland) const {
    if (sigmaRosseland.rank() != 3) {
      OPAC_ERROR("photons::MeanSOpacity: scattering tables must be rank 3");
    }
    if (sigmaRosseland.dim(1) <= 0) {
      OPAC_ERROR("photons::MeanSOpacity: ngroups must be positive");
    }
  }

  template <typename GroupBoundsIndexer>
  void LoadScatteringTables_(const DataBox &sigmaRosseland,
                             const GroupBoundsIndexer &group_bounds) {
    ValidateScatteringTable_(sigmaRosseland);
    ngroups_ = sigmaRosseland.dim(1);
    ValidateGroupBounds(group_bounds, ngroups_);
    SetGroupBounds(groupBounds_, group_bounds, ngroups_);
    lsigmaRosseland_.copyMetadata(sigmaRosseland);
    for (int i = 0; i < sigmaRosseland.size(); ++i) {
      lsigmaRosseland_(i) = ToLog(sigmaRosseland(i));
    }
  }

  void ExportScatteringTables_(DataBox &sigmaRosseland) const {
    sigmaRosseland.copyMetadata(lsigmaRosseland_);
    for (int i = 0; i < lsigmaRosseland_.size(); ++i) {
      sigmaRosseland(i) = FromLog(lsigmaRosseland_(i));
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
    lsigmaRosseland_.resize(NRho, NT, ngroups_);
    lsigmaRosseland_.setRange(1, lTMin, lTMax, NT);
    lsigmaRosseland_.setRange(2, lRhoMin, lRhoMax, NRho);

    PlanckDistribution<PC> dist;
    std::vector<Real> rosselandDenom(ngroups_, 0.);

    for (int iT = 0; iT < NT; ++iT) {
      const Real lT = lsigmaRosseland_.range(1).x(iT);
      const Real T = FromLog(lT);

      for (int group = 0; group < ngroups_; ++group) {
        Real dBdTaccum = 0.;
        const Real nuMin = GroupBoundAt(group_bounds, group);
        const Real nuMax = GroupBoundAt(group_bounds, group + 1);
        ForEachGroupFrequencySample<PC>(
            T, nuMin, nuMax, NNuPerGroup, [&](const Real nu, const Real dnu) {
              Real B = 0.;
              Real dBdT = 0.;
              ThermalWeightsAtNu<PC>(dist, T, nu, B, dBdT);
              dBdTaccum += dBdT * dnu;
            });

        rosselandDenom[group] = dBdTaccum;
      }

      for (int iRho = 0; iRho < NRho; ++iRho) {
        const Real lRho = lsigmaRosseland_.range(2).x(iRho);
        const Real rho = FromLog(lRho);

        for (int group = 0; group < ngroups_; ++group) {
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

                if (sigma > singularity_opac::robust::SMALL()) {
                  sigmaRosselandNum +=
                      singularity_opac::robust::ratio(rho, sigma) * dBdT * dnu;
                }
              });

          const Real sigmaRosseland =
              (rosselandDenom[group] > singularity_opac::robust::SMALL() &&
               sigmaRosselandNum > singularity_opac::robust::SMALL())
                  ? singularity_opac::robust::ratio(rosselandDenom[group],
                                                    sigmaRosselandNum)
                  : 0.;

          lsigmaRosseland_(iRho, iT, group) = ToLog(sigmaRosseland);
          if (std::isnan(lsigmaRosseland_(iRho, iT, group))) {
            OPAC_ERROR("photons::MeanSOpacity: NAN in opacity evaluations");
          }
        }
      }
    }
  }

  DataBox lsigmaRosseland_;
  DataBox groupBounds_;
  int ngroups_ = 0;
};

} // namespace impl

using MeanSOpacityBase = impl::MeanSOpacity<PhysicalConstantsCGS>;

} // namespace photons
} // namespace singularity

#endif // SINGULARITY_OPAC_PHOTONS_MEAN_S_OPACITY_PHOTONS_
