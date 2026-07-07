// ======================================================================
// © 2026. Triad National Security, LLC. All rights reserved.  This
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
#include <singularity-opac/photons/mean_photon_variant.hpp>
#include <singularity-opac/photons/non_cgs_photons.hpp>
#include <singularity-opac/photons/thermal_distributions_photons.hpp>
#include <spiner/databox.hpp>

namespace singularity {
namespace photons {
namespace impl {

#define EPS (10.0 * std::numeric_limits<Real>::min())

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
                    const int NRho, const Real lTMin, const Real lTMax,
                    const int NT, const GroupBoundsIndexer &group_bounds,
                    const int ngroups, const int NNuPerGroup = 64,
                    Real *lambda = nullptr) {
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
    ExportGroupBounds_(groupBounds);

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

  PORTABLE_INLINE_FUNCTION
  Real PlanckMeanAbsorptionCoefficient(const Real rho, const Real temp) const {
    PORTABLE_REQUIRE(ngroups_ == 1,
      "PlanckMeanAbsorptionCoefficient only valid for ngroups==1. "
      "Use PlanckGroupAbsorptionCoefficient(rho, temp, group) for multigroup.");
    return PlanckGroupAbsorptionCoefficient(rho, temp, 0);
  }

  PORTABLE_INLINE_FUNCTION
  Real RosselandMeanAbsorptionCoefficient(const Real rho,
                                          const Real temp) const {
    PORTABLE_REQUIRE(ngroups_ == 1,
      "RosselandMeanAbsorptionCoefficient only valid for ngroups==1. "
      "Use RosselandGroupAbsorptionCoefficient(rho, temp, group) for multigroup.");
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
    if (!(nu >= GroupBoundAt_(groupBounds_, 0) &&
          nu <= GroupBoundAt_(groupBounds_, ngroups_))) {
      OPAC_ERROR(
          "photons::MeanOpacity: frequency is outside group bounds");
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
    const Real lRho = toLog_(rho);
    const Real lT = toLog_(temp);
    return rho * fromLog_(lkappa.interpToReal(lRho, lT, group));
  }

  template <typename GroupBoundsIndexer>
  PORTABLE_INLINE_FUNCTION Real
  GroupBoundAt_(const GroupBoundsIndexer &group_bounds, const int group) const {
    return group_bounds[group];
  }

  PORTABLE_INLINE_FUNCTION
  Real GroupBoundAt_(const DataBox &group_bounds, const int group) const {
    return group_bounds(group);
  }

  template <typename GroupBoundsIndexer>
  void ValidateGroupBounds_(const GroupBoundsIndexer &group_bounds,
                            const int ngroups) const {
    if (ngroups <= 0) {
      OPAC_ERROR("photons::MeanOpacity: ngroups must be positive");
    }
    for (int group = 0; group <= ngroups; ++group) {
      const Real bound = GroupBoundAt_(group_bounds, group);
      if (std::isnan(bound)) {
        OPAC_ERROR("photons::MeanOpacity: group bounds must be finite "
                   "or IEEE +infinity");
      }
      if (std::isinf(bound) && bound < 0.) {
        OPAC_ERROR(
            "photons::MeanOpacity: group bounds may not be -infinity");
      }
      if (group == 0) {
        if (!(bound >= 0.)) {
          OPAC_ERROR("photons::MeanOpacity: first group bound must be "
                     "nonnegative");
        }
      } else if (!(bound > GroupBoundAt_(group_bounds, group - 1))) {
        OPAC_ERROR("photons::MeanOpacity: group bounds must be strictly "
                   "increasing");
      }
      if (!std::isfinite(bound) && group != ngroups) {
        OPAC_ERROR(
            "photons::MeanOpacity: only the final group bound may be "
            "IEEE +infinity");
      }
    }
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
    ValidateGroupBounds_(group_bounds, ngroups_);
    SetGroupBounds_(group_bounds, ngroups_);
    lkappaPlanck_.copyMetadata(kappaPlanck);
    lkappaRosseland_.copyMetadata(kappaRosseland);
    for (int i = 0; i < kappaPlanck.size(); ++i) {
      lkappaPlanck_(i) = toLog_(kappaPlanck(i));
      lkappaRosseland_(i) = toLog_(kappaRosseland(i));
    }
  }

  void ExportOpacityTables_(DataBox &kappaPlanck,
                            DataBox &kappaRosseland) const {
    kappaPlanck.copyMetadata(lkappaPlanck_);
    kappaRosseland.copyMetadata(lkappaRosseland_);
    for (int i = 0; i < lkappaPlanck_.size(); ++i) {
      kappaPlanck(i) = fromLog_(lkappaPlanck_(i));
      kappaRosseland(i) = fromLog_(lkappaRosseland_(i));
    }
  }

  void ExportGroupBounds_(DataBox &groupBounds) const {
    groupBounds.resize(ngroups_ + 1);
    for (int group = 0; group <= ngroups_; ++group) {
      groupBounds(group) = groupBounds_(group);
    }
  }

  template <typename GroupBoundsIndexer>
  void SetGroupBounds_(const GroupBoundsIndexer &group_bounds,
                       const int ngroups) {
    groupBounds_.resize(ngroups + 1);
    for (int group = 0; group <= ngroups; ++group) {
      groupBounds_(group) = GroupBoundAt_(group_bounds, group);
    }
  }

  template <typename SampleOp>
  void ForEachGroupFrequencySample_(const Real temp, const Real nuMin,
                                    const Real nuMax, const int NNuPerGroup,
                                    SampleOp &&sample_op) const {
    // For [0, ∞) or very wide ranges, use thermal-aware sampling
    // For reasonable finite ranges, integrate over the full group bounds
    const Real nu_thermal_min = 1.e-3 * PC::kb * temp / PC::h;
    const Real nu_thermal_max = 1.e3 * PC::kb * temp / PC::h;

    // Determine if we need special handling
    const bool is_lower_extreme = (nuMin == 0.) || (nuMin < 0.1 * nu_thermal_min);
    const bool is_upper_extreme = !std::isfinite(nuMax) || (nuMax > 10. * nu_thermal_max);

    // Set integration bounds, but ensure they're valid
    Real nu_sample_min = is_lower_extreme ? nu_thermal_min : nuMin;
    Real nu_sample_max = is_upper_extreme ? nu_thermal_max : nuMax;

    // If thermal-aware bounds are invalid, use a small but valid range within group bounds
    if (nu_sample_min >= nu_sample_max) {
      if (std::isfinite(nuMax) && nuMax > 0.) {
        // Group is [0 or small, nuMax]: sample near nuMax
        nu_sample_min = 0.5 * nuMax;
        nu_sample_max = nuMax;
      } else {
        // Group extends to infinity: sample around thermal peak
        nu_sample_min = 0.1 * nu_thermal_max;
        nu_sample_max = nu_thermal_max;
      }
    }

    // Use logarithmic spacing with midpoint rule
    const Real lNuMin = toLog_(nu_sample_min);
    const Real lNuMax = toLog_(nu_sample_max);
    const Real dlnu = (lNuMax - lNuMin) / NNuPerGroup;
    for (int inu = 0; inu < NNuPerGroup; ++inu) {
      const Real lnu = lNuMin + (inu + 0.5) * dlnu;
      const Real nu = fromLog_(lnu);
      sample_op(nu, nu * dlnu);
    }
  }

  void ThermalWeightsAtNu_(const PlanckDistribution<PC> &dist, const Real temp,
                           const Real nu, Real &B, Real &dBdT) const {
    const Real x = PC::h * nu / (PC::kb * temp);
    if (x < 80.) {
      B = dist.ThermalDistributionOfTNu(temp, nu);
      dBdT = dist.DThermalDistributionOfTNuDT(temp, nu);
      return;
    }

    const Real expMinusX = std::exp(-x);
    B = (2. * PC::h * nu * nu * nu / (PC::c * PC::c)) * expMinusX;
    dBdT = 2. * PC::h * PC::h * nu * nu * nu * nu * expMinusX /
           (temp * temp * PC::c * PC::c * PC::kb);
  }

  template <typename Opacity, typename GroupBoundsIndexer>
  void MeanOpacityImpl_(const Opacity &opac, const Real lRhoMin,
                              const Real lRhoMax, const int NRho,
                              const Real lTMin, const Real lTMax, const int NT,
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
    ValidateGroupBounds_(group_bounds, ngroups);

    ngroups_ = ngroups;
    SetGroupBounds_(group_bounds, ngroups_);
    lkappaPlanck_.resize(NRho, NT, ngroups_);
    lkappaPlanck_.setRange(1, lTMin, lTMax, NT);
    lkappaPlanck_.setRange(2, lRhoMin, lRhoMax, NRho);
    lkappaRosseland_.copyMetadata(lkappaPlanck_);

    PlanckDistribution<PC> dist;
    std::vector<Real> planckDenom(ngroups_, 0.);
    std::vector<Real> rosselandDenom(ngroups_, 0.);

    for (int iT = 0; iT < NT; ++iT) {
      const Real lT = lkappaPlanck_.range(1).x(iT);
      const Real T = fromLog_(lT);

      for (int group = 0; group < ngroups_; ++group) {
        Real Baccum = 0.;
        Real dBdTaccum = 0.;
        const Real nuMin = GroupBoundAt_(group_bounds, group);
        const Real nuMax = GroupBoundAt_(group_bounds, group + 1);
        ForEachGroupFrequencySample_(
            T, nuMin, nuMax, NNuPerGroup, [&](const Real nu, const Real dnu) {
              Real B = 0.;
              Real dBdT = 0.;
              ThermalWeightsAtNu_(dist, T, nu, B, dBdT);
              Baccum += B * dnu;
              dBdTaccum += dBdT * dnu;
            });

        planckDenom[group] = Baccum;
        rosselandDenom[group] = dBdTaccum;
      }

      for (int iRho = 0; iRho < NRho; ++iRho) {
        const Real lRho = lkappaPlanck_.range(2).x(iRho);
        const Real rho = fromLog_(lRho);

        for (int group = 0; group < ngroups_; ++group) {
          Real kappaPlanckNum = 0.;
          Real kappaRosselandNum = 0.;
          const Real nuMin = GroupBoundAt_(group_bounds, group);
          const Real nuMax = GroupBoundAt_(group_bounds, group + 1);
          ForEachGroupFrequencySample_(
              T, nuMin, nuMax, NNuPerGroup, [&](const Real nu, const Real dnu) {
                const Real alpha =
                    opac.AbsorptionCoefficient(rho, T, nu, lambda);
                Real B = 0.;
                Real dBdT = 0.;
                ThermalWeightsAtNu_(dist, T, nu, B, dBdT);
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

          lkappaPlanck_(iRho, iT, group) = toLog_(kappaPlanck);
          lkappaRosseland_(iRho, iT, group) = toLog_(kappaRosseland);
          if (std::isnan(lkappaPlanck_(iRho, iT, group)) ||
              std::isnan(lkappaRosseland_(iRho, iT, group))) {
            OPAC_ERROR(
                "photons::MeanOpacity: NAN in opacity evaluations");
          }
        }
      }
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real toLog_(const Real x) const { return std::log10(std::abs(x) + EPS); }

  PORTABLE_INLINE_FUNCTION
  Real fromLog_(const Real lx) const { return std::pow(10., lx); }

  PORTABLE_INLINE_FUNCTION
  int GroupOfNuImpl_(const Real nu) const {
    if (nu == GroupBoundAt_(groupBounds_, ngroups_)) {
      return ngroups_ - 1;
    }
    int lower = 0;
    int upper = ngroups_;
    while (upper - lower > 1) {
      const int middle = (lower + upper) / 2;
      if (nu < GroupBoundAt_(groupBounds_, middle)) {
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

#undef EPS

} // namespace impl

using MeanOpacityBase = impl::MeanOpacity<PhysicalConstantsCGS>;
using MeanOpacity =
    impl::MeanVariant<MeanOpacityBase,
                            MeanNonCGSUnits<MeanOpacityBase>>;

} // namespace photons
} // namespace singularity

#endif // SINGULARITY_OPAC_PHOTONS_MEAN_OPACITY_PHOTONS_
