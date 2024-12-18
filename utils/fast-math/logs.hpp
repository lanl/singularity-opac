//======================================================================
// approximate log10, approximately 5x faster than std::log10
// Author: Jonah Miller (jonahm@lanl.gov)

// © 2021. Triad National Security, LLC. All rights reserved.  This
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

#ifndef _FAST_MATH_LOGS_
#define _FAST_MATH_LOGS_

#include <cmath>
#include <ports-of-call/portability.hpp>

// TODO(JMM): This has got to go. Replace with sane math.

#ifdef SINGULARITY_USE_FMATH
// herumi-fmath does not work on device
// On CPUS it provides another 10% or so speedup
// for negligible cost in accuracy
#define BD_USE_FMATH
#if defined(PORTABILITY_STRATEGY_KOKKOS)
  #undef BD_USE_FMATH
#else
  #include <fmath.hpp>
#endif
#endif // SINGULARITY_USE_FMATH

namespace BDMath {

PORTABLE_INLINE_FUNCTION
float log10(const float x) {
  // const double ILOG10 = 1./std::log(10.0);
  constexpr double ILOG10 = 0.43429448190325176;
#ifdef BD_USE_FMATH
  return fmath::log(x) * ILOG10;
#else
  return logf(x) * ILOG10;
#endif
}

} // namespace BDMath

#undef BD_USE_FMATH
#endif //  _FAST_MATH_LOGS_
