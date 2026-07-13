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

#ifndef SINGULARITY_OPAC_PHOTONS_MEAN_PHOTON_UTILS_
#define SINGULARITY_OPAC_PHOTONS_MEAN_PHOTON_UTILS_

// This file was made in part with generative AI.

#include <cmath>
#include <limits>

#include <ports-of-call/portability.hpp>
#include <singularity-opac/base/opac_error.hpp>
#include <singularity-opac/photons/mean_photon_types.hpp>
#include <singularity-opac/photons/thermal_distributions_photons.hpp>
#include <spiner/databox.hpp>

namespace singularity {
namespace photons {
namespace impl {

using MeanUtilsDataBox = Spiner::DataBox<Real>;

// Log/anti-log transforms used to store and interpolate opacities. A small
// floor keeps toLog well-defined at zero.
PORTABLE_INLINE_FUNCTION
Real ToLog(const Real x) {
  constexpr Real eps = 10.0 * std::numeric_limits<Real>::min();
  return std::log10(std::abs(x) + eps);
}

PORTABLE_INLINE_FUNCTION
Real FromLog(const Real lx) { return std::pow(10., lx); }

// Group-bound access, supporting both a generic indexer (used at construction
// time) and a stored DataBox (used after the bounds are cached).
template <typename GroupBoundsIndexer>
PORTABLE_INLINE_FUNCTION Real
GroupBoundAt(const GroupBoundsIndexer &group_bounds, const int group) {
  return group_bounds[group];
}

PORTABLE_INLINE_FUNCTION
Real GroupBoundAt(const MeanUtilsDataBox &group_bounds, const int group) {
  return group_bounds(group);
}

// Copy user-provided group bounds into the class-owned DataBox.
template <typename GroupBoundsIndexer>
void SetGroupBounds(MeanUtilsDataBox &groupBounds,
                    const GroupBoundsIndexer &group_bounds, const int ngroups) {
  groupBounds.resize(ngroups + 1);
  for (int group = 0; group <= ngroups; ++group) {
    groupBounds(group) = GroupBoundAt(group_bounds, group);
  }
}

// Copy the class-owned group bounds back out (e.g. for HDF5 export).
inline void ExportGroupBounds(MeanUtilsDataBox &groupBounds,
                              const MeanUtilsDataBox &storedBounds,
                              const int ngroups) {
  groupBounds.resize(ngroups + 1);
  for (int group = 0; group <= ngroups; ++group) {
    groupBounds(group) = storedBounds(group);
  }
}

// Validate half-open group bounds [nu_g, nu_{g+1}): strictly increasing,
// nonnegative, with only the final bound permitted to be IEEE +infinity.
template <typename GroupBoundsIndexer>
void ValidateGroupBounds(const GroupBoundsIndexer &group_bounds,
                         const int ngroups) {
  if (ngroups <= 0) {
    OPAC_ERROR("photons multigroup: ngroups must be positive");
  }
  for (int group = 0; group <= ngroups; ++group) {
    const Real bound = GroupBoundAt(group_bounds, group);
    if (std::isnan(bound)) {
      OPAC_ERROR("photons multigroup: group bounds must be finite "
                 "or IEEE +infinity");
    }
    if (std::isinf(bound) && bound < 0.) {
      OPAC_ERROR("photons multigroup: group bounds may not be -infinity");
    }
    if (group == 0) {
      if (!(bound >= 0.)) {
        OPAC_ERROR("photons multigroup: first group bound must be "
                   "nonnegative");
      }
    } else if (!(bound > GroupBoundAt(group_bounds, group - 1))) {
      OPAC_ERROR("photons multigroup: group bounds must be strictly "
                 "increasing");
    }
    if (!std::isfinite(bound) && group != ngroups) {
      OPAC_ERROR("photons multigroup: only the final group bound may be "
                 "IEEE +infinity");
    }
  }
}

// Sample a group's frequency range on a logarithmic, midpoint-rule grid,
// invoking sample_op(nu, weight) for each of NNuPerGroup points. Extreme
// bounds (nu=0, +infinity, or far from the thermal peak) are clamped to a
// thermal-aware window so the integral stays well-conditioned.
template <typename PC, typename SampleOp>
void ForEachGroupFrequencySample(const Real temp, const Real nuMin,
                                 const Real nuMax, const int NNuPerGroup,
                                 SampleOp &&sample_op) {
  // For [0, ∞) or very wide ranges, use thermal-aware sampling
  // For reasonable finite ranges, integrate over the full group bounds
  const Real nu_thermal_min = 1.e-3 * PC::kb * temp / PC::h;
  const Real nu_thermal_max = 1.e3 * PC::kb * temp / PC::h;

  // Determine if we need special handling
  const bool is_lower_extreme = (nuMin == 0.) || (nuMin < 0.1 * nu_thermal_min);
  const bool is_upper_extreme =
      !std::isfinite(nuMax) || (nuMax > 10. * nu_thermal_max);

  // Set integration bounds, but ensure they're valid
  Real nu_sample_min = is_lower_extreme ? nu_thermal_min : nuMin;
  Real nu_sample_max = is_upper_extreme ? nu_thermal_max : nuMax;

  // If thermal-aware bounds are invalid, use a small but valid range within
  // group bounds
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
  const Real lNuMin = ToLog(nu_sample_min);
  const Real lNuMax = ToLog(nu_sample_max);
  const Real dlnu = (lNuMax - lNuMin) / NNuPerGroup;
  for (int inu = 0; inu < NNuPerGroup; ++inu) {
    const Real lnu = lNuMin + (inu + 0.5) * dlnu;
    const Real nu = FromLog(lnu);
    sample_op(nu, nu * dlnu);
  }
}

// Evaluate the Planck function B_nu and its temperature derivative dB_nu/dT,
// switching to the Wien closed form deep in the tail (see wien_tail_x).
template <typename PC>
void ThermalWeightsAtNu(const PlanckDistribution<PC> &dist, const Real temp,
                        const Real nu, Real &B, Real &dBdT) {
  const Real x = PC::h * nu / (PC::kb * temp);
  if (x < wien_tail_x) {
    B = dist.ThermalDistributionOfTNu(temp, nu);
    dBdT = dist.DThermalDistributionOfTNuDT(temp, nu);
    return;
  }

  const Real expMinusX = std::exp(-x);
  B = (2. * PC::h * nu * nu * nu / (PC::c * PC::c)) * expMinusX;
  dBdT = 2. * PC::h * PC::h * nu * nu * nu * nu * expMinusX /
         (temp * temp * PC::c * PC::c * PC::kb);
}

} // namespace impl
} // namespace photons
} // namespace singularity

#endif // SINGULARITY_OPAC_PHOTONS_MEAN_PHOTON_UTILS_
