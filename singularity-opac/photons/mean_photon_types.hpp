// ======================================================================
// © 2025-2026. Triad National Security, LLC. All rights reserved.  This
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

#ifndef SINGULARITY_OPAC_PHOTONS_MEAN_PHOTON_TYPES_
#define SINGULARITY_OPAC_PHOTONS_MEAN_PHOTON_TYPES_

// This file was made in part with generative AI.

#include <ports-of-call/portability.hpp>

namespace singularity {
namespace photons {

// mean-opacity mode
enum OpacityAveraging { Rosseland = 0, Planck = 1 };

// Crossover into the Wien tail of the Planck function, in units of the
// dimensionless x = h * nu / (kb * temp). For x above this value the exact
// Planck form 1 / expm1(x) and the Wien approximation exp(-x) agree to well
// below machine precision (exp(-80) ~ 2e-35), so we switch to the cheaper,
// monotone Wien closed form. This is a regime switch, not an overflow guard.
constexpr Real wien_tail_x = 80.;

} // namespace photons
} // namespace singularity

#endif // SINGULARITY_OPAC_PHOTONS_MEAN_PHOTON_TYPES_
