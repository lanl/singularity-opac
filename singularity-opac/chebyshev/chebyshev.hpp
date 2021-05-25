// ======================================================================
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

#ifndef SINGULARITY_OPAC_CHEBYSHEV_CHEBYSHEV_
#define SINGULARITY_OPAC_CHEBYSHEV_CHEBYSHEV_

#include <cassert>
#include <cmath>
#include <iostream>

#include <ports-of-call/portability.hpp>

#include <singularity-opac/chebyshev/vandermonde.hpp>

// Routines for Chebyshev interpolation and integration.  For more
// details, see Heath, 1997 or Boydm 1941.

namespace singularity {
namespace chebyshev {

template <typename Matrix, typename VecIn, typename VecOut>
PORTABLE_INLINE_FUNCTION void MatMultiply(const Matrix &a, const VecIn &x,
                                          VecOut &y, int npoints) {
  for (int i = 0; i < npoints; ++i) {
    y[i] = 0;
    for (int j = 0; j < npoints; ++j) {
      y[i] += a[i][j] * x[j];
    }
  }
}

// Unit cell is [-1,1]
PORTABLE_INLINE_FUNCTION Real ToUnitCell(Real x, Real xmin, Real xmax) {
  assert(xmax > xmin);
  return 2 * (x - xmin) / (xmax - xmin) - 1;
}
PORTABLE_INLINE_FUNCTION
Real FromUnitCell(Real x, Real xmin, Real xmax) {
  assert(xmax > xmin);
  return (x + 1) * (xmax - xmin) / 2 + xmin;
}

template <typename Indexer>
PORTABLE_INLINE_FUNCTION void GetPoints(const Real xmin, const Real xmax,
                                        const int npoints, Indexer &points) {
  for (int k = 1; k <= npoints; ++k) {
    Real tk = -std::cos((2 * k - 1) * M_PI / (2 * npoints));
    points[k - 1] = FromUnitCell(tk, xmin, xmax);
  }
}

// TODO(JMM): Should this be templated? Or Recursive? What's faster?
// Templates, iteration, or recursion?
// Recursion relation for Chebyshev polynomials of the first kind:
// T_0(x) = 1
// T_1(x) = x
// T_i(x) = 2*x*T_{i-1}(x) - T_{i-2}(x)
PORTABLE_INLINE_FUNCTION
Real T(const int i, const Real x) {
  Real Ti[3];
  Ti[0] = 1;
  Ti[1] = x;
  Ti[2] = 2 * x * Ti[1] - Ti[0];
  if (i < 3) return Ti[i];
  for (int ii = 3; ii <= i; ++ii) {
    Ti[0] = Ti[1];
    Ti[1] = Ti[2];
    Ti[2] = 2 * x * Ti[1] - Ti[0];
  }
  return Ti[2];
}

// Chebyshev interpolation from coefficients
template <typename Indexer>
PORTABLE_INLINE_FUNCTION Real InterpFromCoeffs(Real x, const Real xmin,
                                               const Real xmax,
                                               const Indexer &coeffs,
                                               const int ncoeffs) {
  Real out = 0;
  x = ToUnitCell(x, xmin, xmax);
  for (int i = 0; i < ncoeffs; ++i) {
    out += coeffs[i] * T(i, x);
  }
  return out;
}

} // namespace chebyshev
} // namespace singularity

#endif // SINGULARITY_OPAC_CHEBYSHEV_CHEBYSHEV_
