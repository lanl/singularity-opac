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
    LoadOpacityTables_(kappaPlanck, kappaRosseland, group_bounds);
  }

#ifdef SPINER_USE_HDF
  MeanOpacity(const std::string &filename) {
    // HDF-backed multigroup tables are expected to provide an ngroups + 1
    // "group bounds" dataset. If the last group is [nu_max, infinity), then
    // the final stored bound must be IEEE +infinity.
    DataBox kappaPlanck;
    DataBox kappaRosseland;
    DataBox groupBounds;
    herr_t status = H5_SUCCESS;
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    status +=
        kappaPlanck.loadHDF(file, SP5::MultigroupOpac::PlanckGroupOpacity);
    status += kappaRosseland.loadHDF(
        file, SP5::MultigroupOpac::RosselandGroupOpacity);
    status += groupBounds.loadHDF(file, SP5::MultigroupOpac::GroupBounds);
    status += H5Fclose(file);

    if (status != H5_SUCCESS) {
      OPAC_ERROR("photons::MeanOpacity: HDF5 error\n");
    }

    LoadOpacityTables_(kappaPlanck, kappaRosseland, groupBounds);
    groupBounds.finalize();
    kappaPlanck.finalize();
    kappaRosseland.finalize();
  }

  void Save(const std::string &filename) const {
    DataBox kappaPlanck;
    DataBox kappaRosseland;
    DataBox groupBounds;
    ExportOpacityTables_(kappaPlanck, kappaRosseland);
    ExportGroupBounds(groupBounds, groupBounds_, ngroups_);

    herr_t status = H5_SUCCESS;
    hid_t file =
        H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    status +=
        kappaPlanck.saveHDF(file, SP5::MultigroupOpac::PlanckGroupOpacity);
    status += kappaRosseland.saveHDF(
        file, SP5::MultigroupOpac::RosselandGroupOpacity);
    status += groupBounds.saveHDF(file, SP5::MultigroupOpac::GroupBounds);
    status += H5Fclose(file);

    kappaPlanck.finalize();
    kappaRosseland.finalize();
    groupBounds.finalize();

    if (status != H5_SUCCESS) {
      OPAC_ERROR("photons::MeanOpacity: HDF5 error\n");
    }
  }
#endif

  PORTABLE_INLINE_FUNCTION
  void PrintParams() const {
    printf("Photon multigroup opacity. ngroups = %d\n", ngroups_);
  }

  MeanOpacity GetOnDevice() {
    MeanOpacity other;
    other.lkappaPlanck_ = Spiner::getOnDeviceDataBox(lkappaPlanck_);
    other.lkappaRosseland_ = Spiner::getOnDeviceDataBox(lkappaRosseland_);
    other.groupBounds_ = Spiner::getOnDeviceDataBox(groupBounds_);
    other.ngroups_ = ngroups_;
    return other;
  }

  void Finalize() {
    lkappaPlanck_.finalize();
    lkappaRosseland_.finalize();
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
    return GroupAbsorptionCoefficient_(lkappaPlanck_, rho, temp, group);
  }

  PORTABLE_INLINE_FUNCTION
  Real RosselandGroupAbsorptionCoefficient(const Real rho, const Real temp,
                                           const int group) const {
    return GroupAbsorptionCoefficient_(lkappaRosseland_, rho, temp, group);
  }

  PORTABLE_INLINE_FUNCTION
  Real AbsorptionCoefficient(const Real rho, const Real temp, const int group,
                             const int gmode = Rosseland) const {
    return (gmode == Planck)
               ? PlanckGroupAbsorptionCoefficient(rho, temp, group)
               : RosselandGroupAbsorptionCoefficient(rho, temp, group);
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
    return GroupOfNuImpl_(nu);
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
  PORTABLE_INLINE_FUNCTION
  Real GroupAbsorptionCoefficient_(const DataBox &lkappa, const Real rho,
                                   const Real temp, const int group) const {
    const Real lRho = ToLog(rho);
    const Real lT = ToLog(temp);
    return rho * FromLog(lkappa.interpToReal(lRho, lT, group));
  }

  void ValidateOpacityTables_(const DataBox &kappaPlanck,
                              const DataBox &kappaRosseland) const {
    if (kappaPlanck.rank() != 3 || kappaRosseland.rank() != 3) {
      OPAC_ERROR("photons::MeanOpacity: opacity tables must be rank 3");
    }
    for (int dim = 1; dim <= 3; ++dim) {
      if (kappaPlanck.dim(dim) != kappaRosseland.dim(dim)) {
        OPAC_ERROR("photons::MeanOpacity: table dimensions do not match");
      }
    }
    if (kappaPlanck.dim(1) <= 0) {
      OPAC_ERROR("photons::MeanOpacity: ngroups must be positive");
    }
    if (kappaPlanck.range(1) != kappaRosseland.range(1) ||
        kappaPlanck.range(2) != kappaRosseland.range(2)) {
      OPAC_ERROR("photons::MeanOpacity: table ranges do not match");
    }
  }

  template <typename GroupBoundsIndexer>
  void LoadOpacityTables_(const DataBox &kappaPlanck,
                          const DataBox &kappaRosseland,
                          const GroupBoundsIndexer &group_bounds) {
    ValidateOpacityTables_(kappaPlanck, kappaRosseland);
    ngroups_ = kappaPlanck.dim(1);
    ValidateGroupBounds(group_bounds, ngroups_);
    SetGroupBounds(groupBounds_, group_bounds, ngroups_);
    lkappaPlanck_.copyMetadata(kappaPlanck);
    lkappaRosseland_.copyMetadata(kappaRosseland);
    for (int i = 0; i < kappaPlanck.size(); ++i) {
      lkappaPlanck_(i) = ToLog(kappaPlanck(i));
      lkappaRosseland_(i) = ToLog(kappaRosseland(i));
    }
  }

  void ExportOpacityTables_(DataBox &kappaPlanck,
                            DataBox &kappaRosseland) const {
    kappaPlanck.copyMetadata(lkappaPlanck_);
    kappaRosseland.copyMetadata(lkappaRosseland_);
    for (int i = 0; i < lkappaPlanck_.size(); ++i) {
      kappaPlanck(i) = FromLog(lkappaPlanck_(i));
      kappaRosseland(i) = FromLog(lkappaRosseland_(i));
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

  PORTABLE_INLINE_FUNCTION
  int GroupOfNuImpl_(const Real nu) const {
    // Shortcuts for boundary cases
    if (nu <= GroupBoundAt(groupBounds_, 0)) {
      return 0;
    }
    if (nu >= GroupBoundAt(groupBounds_, ngroups_)) {
      return ngroups_ - 1;
    }
    // Binary search to find group index containing nu
    int lower = 0;
    int upper = ngroups_;
    while (upper - lower > 1) {
      const int middle = (lower + upper) / 2;
      if (nu < GroupBoundAt(groupBounds_, middle)) {
        upper = middle;
      } else {
        lower = middle;
      }
    }
    return lower;
  }

  DataBox lkappaPlanck_;
  DataBox lkappaRosseland_;
  DataBox groupBounds_;
  int ngroups_ = 0;
};

} // namespace impl

using MeanOpacityBase = impl::MeanOpacity<PhysicalConstantsCGS>;
using MeanOpacity =
    impl::MeanVariant<MeanOpacityBase, MeanNonCGSUnits<MeanOpacityBase>>;

} // namespace photons
} // namespace singularity

#endif // SINGULARITY_OPAC_PHOTONS_MEAN_OPACITY_PHOTONS_
