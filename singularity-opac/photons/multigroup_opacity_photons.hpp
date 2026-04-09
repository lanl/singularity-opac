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

#ifndef SINGULARITY_OPAC_PHOTONS_MULTIGROUP_OPACITY_PHOTONS_
#define SINGULARITY_OPAC_PHOTONS_MULTIGROUP_OPACITY_PHOTONS_

// This file was made in part with generative AI.

#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include <ports-of-call/portability.hpp>
#include <singularity-opac/base/opac_error.hpp>
#include <singularity-opac/base/robust_utils.hpp>
#include <singularity-opac/base/sp5.hpp>
#include <singularity-opac/constants/constants.hpp>
#include <singularity-opac/photons/mean_photon_types.hpp>
#include <singularity-opac/photons/multigroup_photon_variant.hpp>
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
class MultigroupOpacity {
 public:
  using PC = pc;
  using DataBox = Spiner::DataBox<Real>;

  MultigroupOpacity() = default;

  template <typename Opacity, typename GroupBoundsIndexer>
  MultigroupOpacity(const Opacity &opac, const Real lRhoMin, const Real lRhoMax,
                    const int NRho, const Real lTMin, const Real lTMax,
                    const int NT, const GroupBoundsIndexer &group_bounds,
                    const int ngroups, const int NNuPerGroup = 64,
                    Real *lambda = nullptr) {
    MultigroupOpacityImpl_(opac, lRhoMin, lRhoMax, NRho, lTMin, lTMax, NT,
                           group_bounds, ngroups, NNuPerGroup, lambda);
  }

  template <typename Opacity>
  MultigroupOpacity(const Opacity &opac, const Real lRhoMin, const Real lRhoMax,
                    const int NRho, const Real lTMin, const Real lTMax,
                    const int NT, const Real nuMin, const Real nuMax,
                    const int NLogGroups, const int NNuPerGroup = 64,
                    Real *lambda = nullptr) {
    auto group_bounds =
        LogSpacedGroupBoundsWithTailGroups_(nuMin, nuMax, NLogGroups);
    MultigroupOpacityImpl_(opac, lRhoMin, lRhoMax, NRho, lTMin, lTMax, NT,
                           group_bounds, TotalGroupsWithTails_(NLogGroups),
                           NNuPerGroup, lambda);
  }

  template <typename GroupBoundsIndexer>
  MultigroupOpacity(const DataBox &kappaPlanck, const DataBox &kappaRosseland,
                    const GroupBoundsIndexer &group_bounds) {
    // Table-backed multigroup opacities always carry explicit group bounds.
    // To represent [nu_max, infinity), the final bound must literally be
    // IEEE +infinity, not a large finite proxy value.
    LoadOpacityTables_(kappaPlanck, kappaRosseland, group_bounds);
  }

#ifdef SPINER_USE_HDF
  MultigroupOpacity(const std::string &filename) {
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
      OPAC_ERROR("photons::MultigroupOpacity: HDF5 error\n");
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
      OPAC_ERROR("photons::MultigroupOpacity: HDF5 error\n");
    }
  }
#endif

  PORTABLE_INLINE_FUNCTION
  void PrintParams() const {
    printf("Photon multigroup opacity. ngroups = %d\n", ngroups_);
  }

  MultigroupOpacity GetOnDevice() {
    MultigroupOpacity other;
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
  int GroupOfNu(const Real nu) const {
    if (!(nu >= GroupBoundAt_(groupBounds_, 0) &&
          nu <= GroupBoundAt_(groupBounds_, ngroups_))) {
      OPAC_ERROR(
          "photons::MultigroupOpacity: frequency is outside group bounds");
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
      OPAC_ERROR("photons::MultigroupOpacity: ngroups must be positive");
    }
    for (int group = 0; group <= ngroups; ++group) {
      const Real bound = GroupBoundAt_(group_bounds, group);
      if (std::isnan(bound)) {
        OPAC_ERROR("photons::MultigroupOpacity: group bounds must be finite "
                   "or IEEE +infinity");
      }
      if (std::isinf(bound) && bound < 0.) {
        OPAC_ERROR(
            "photons::MultigroupOpacity: group bounds may not be -infinity");
      }
      if (group == 0) {
        if (!(bound >= 0.)) {
          OPAC_ERROR("photons::MultigroupOpacity: first group bound must be "
                     "nonnegative");
        }
      } else if (!(bound > GroupBoundAt_(group_bounds, group - 1))) {
        OPAC_ERROR("photons::MultigroupOpacity: group bounds must be strictly "
                   "increasing");
      }
      if (!std::isfinite(bound) && group != ngroups) {
        OPAC_ERROR(
            "photons::MultigroupOpacity: only the final group bound may be "
            "IEEE +infinity");
      }
    }
  }

  void ValidateLogGroupInputs_(const Real nuMin, const Real nuMax,
                               const int NLogGroups) const {
    if (!(nuMin > 0.)) {
      OPAC_ERROR("photons::MultigroupOpacity: nuMin must be positive");
    }
    if (!(nuMax > nuMin)) {
      OPAC_ERROR("photons::MultigroupOpacity: nuMax must be greater than "
                 "nuMin");
    }
    if (NLogGroups <= 0) {
      OPAC_ERROR("photons::MultigroupOpacity: NLogGroups must be positive");
    }
  }

  int TotalGroupsWithTails_(const int NLogGroups) const {
    return NLogGroups + 2;
  }

  std::vector<Real>
  LogSpacedGroupBoundsWithTailGroups_(const Real nuMin, const Real nuMax,
                                      const int NLogGroups) const {
    ValidateLogGroupInputs_(nuMin, nuMax, NLogGroups);
    const int ngroups = TotalGroupsWithTails_(NLogGroups);
    std::vector<Real> group_bounds(ngroups + 1, 0.);
    group_bounds[0] = 0.;
    group_bounds[1] = nuMin;
    group_bounds[ngroups - 1] = nuMax;
    group_bounds[ngroups] = std::numeric_limits<Real>::infinity();

    const Real lNuMin = std::log10(nuMin);
    const Real lNuMax = std::log10(nuMax);
    for (int group = 1; group < NLogGroups; ++group) {
      const Real frac = static_cast<Real>(group) / NLogGroups;
      group_bounds[group + 1] =
          std::pow(10., lNuMin + frac * (lNuMax - lNuMin));
    }
    return group_bounds;
  }

  void ValidateOpacityTables_(const DataBox &kappaPlanck,
                              const DataBox &kappaRosseland) const {
    if (kappaPlanck.rank() != 3 || kappaRosseland.rank() != 3) {
      OPAC_ERROR("photons::MultigroupOpacity: opacity tables must be rank 3");
    }
    for (int dim = 1; dim <= 3; ++dim) {
      if (kappaPlanck.dim(dim) != kappaRosseland.dim(dim)) {
        OPAC_ERROR("photons::MultigroupOpacity: table dimensions do not match");
      }
    }
    if (kappaPlanck.dim(1) <= 0) {
      OPAC_ERROR("photons::MultigroupOpacity: ngroups must be positive");
    }
    if (kappaPlanck.range(1) != kappaRosseland.range(1) ||
        kappaPlanck.range(2) != kappaRosseland.range(2)) {
      OPAC_ERROR("photons::MultigroupOpacity: table ranges do not match");
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
    if (nuMin == 0.) {
      const Real du = 1. / NNuPerGroup;
      for (int inu = 0; inu < NNuPerGroup; ++inu) {
        const Real u = (inu + 0.5) * du;
        const Real nu = nuMax * u * u;
        const Real dnu = 2. * nuMax * u * du;
        sample_op(nu, dnu);
      }
      return;
    }

    if (!std::isfinite(nuMax)) {
      const Real du = 1. / NNuPerGroup;
      const Real xMin = PC::h * nuMin / (PC::kb * temp);
      const Real nuScale = PC::kb * temp / PC::h;
      for (int inu = 0; inu < NNuPerGroup; ++inu) {
        const Real u = (inu + 0.5) * du;
        const Real oneMinusU = 1. - u;
        const Real x = xMin + u / oneMinusU;
        const Real nu = nuScale * x;
        const Real dnu = nuScale * du / (oneMinusU * oneMinusU);
        sample_op(nu, dnu);
      }
      return;
    }

    const Real lNuMin = toLog_(nuMin);
    const Real lNuMax = toLog_(nuMax);
    const Real dlnu = (lNuMax - lNuMin) / (NNuPerGroup - 1);
    for (int inu = 0; inu < NNuPerGroup; ++inu) {
      const Real weight = (inu == 0 || inu == NNuPerGroup - 1) ? 0.5 : 1.;
      const Real lnu = lNuMin + inu * dlnu;
      const Real nu = fromLog_(lnu);
      sample_op(nu, weight * nu * dlnu);
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
  void MultigroupOpacityImpl_(const Opacity &opac, const Real lRhoMin,
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
      OPAC_ERROR("photons::MultigroupOpacity: NNuPerGroup must be at least 2");
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
              rosselandDenom[group] > singularity_opac::robust::SMALL()
                  ? singularity_opac::robust::ratio(rosselandDenom[group],
                                                    kappaRosselandNum)
                  : 0.;

          lkappaPlanck_(iRho, iT, group) = toLog_(kappaPlanck);
          lkappaRosseland_(iRho, iT, group) = toLog_(kappaRosseland);
          if (std::isnan(lkappaPlanck_(iRho, iT, group)) ||
              std::isnan(lkappaRosseland_(iRho, iT, group))) {
            OPAC_ERROR(
                "photons::MultigroupOpacity: NAN in opacity evaluations");
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

using MultigroupOpacityBase = impl::MultigroupOpacity<PhysicalConstantsCGS>;
using MultigroupOpacity =
    impl::MultigroupVariant<MultigroupOpacityBase,
                            MultigroupNonCGSUnits<MultigroupOpacityBase>>;

} // namespace photons
} // namespace singularity

#endif // SINGULARITY_OPAC_PHOTONS_MULTIGROUP_OPACITY_PHOTONS_
