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

#include <ports-of-call/portability.hpp>

// Routines for Chebyshev interpolation and integration.  For more
// details, see Heath, 1997 or Boydm 1941.

namespace singularity {
namespace chebyshev {

// Inverse Vandermonde matrices computed via Mathematica
constexpr Real Vandermonde3[3][3] = {
    {0.33333333333333, 0.33333333333333, 0.33333333333333},
    {-0.57735026918963, 0, 0.57735026918963},
    {0.33333333333333, -0.66666666666667, 0.33333333333333}};
constexpr Real Vandermonde5[5][5] = {
    {0.20000000000000, 0.20000000000000, 0.20000000000000, 0.20000000000000,
     0.20000000000000},
    {-0.38042260651806, -0.23511410091699, 0, 0.23511410091699,
     0.38042260651806},
    {0.32360679774998, -0.12360679774998, -0.40000000000000, -0.12360679774998,
     0.32360679774998},
    {-0.23511410091699, 0.38042260651806, 0, -0.38042260651806,
     0.23511410091699},
    {0.12360679774998, -0.32360679774998, 0.40000000000000, -0.32360679774998,
     0.12360679774998}};
constexpr Real Vandermonde7[7][7] = {
    {0.14285714285714, 0.14285714285714, 0.14285714285714, 0.14285714285714,
     0.14285714285714, 0.14285714285714, 0.14285714285714},
    {-0.27855083205195, -0.22338042356229, -0.12396678260502, 0,
     0.12396678260502, 0.22338042356229, 0.27855083205195},
    {0.25741967654355, 0.063577409701804, -0.17813994338821, -0.28571428571429,
     -0.17813994338821, 0.063577409701804, 0.25741967654355},
    {-0.22338042356229, 0.12396678260502, 0.27855083205195, 0,
     -0.27855083205195, -0.12396678260502, 0.22338042356229},
    {0.17813994338821, -0.25741967654355, -0.063577409701804, 0.28571428571429,
     -0.063577409701804, -0.25741967654355, 0.17813994338821},
    {-0.12396678260502, 0.27855083205195, -0.22338042356229, 0,
     0.22338042356229, -0.27855083205195, 0.12396678260502},
    {0.063577409701804, -0.17813994338821, 0.25741967654355, -0.28571428571429,
     0.25741967654355, -0.17813994338821, 0.063577409701804}};
constexpr Real Vandermonde9[9][9] = {
    {0.11111111111111, 0.11111111111111, 0.11111111111111, 0.11111111111111,
     0.11111111111111, 0.11111111111111, 0.11111111111111, 0.11111111111111,
     0.11111111111111},
    {-0.21884616733605, -0.19245008972988, -0.14284169104145,
     -0.076004476294593, 0, 0.076004476294593, 0.14284169104145,
     0.19245008972988, 0.21884616733605},
    {0.20882058239687, 0.11111111111111, -0.038588483925985, -0.17023209847088,
     -0.22222222222222, -0.17023209847088, -0.038588483925985, 0.11111111111111,
     0.20882058239687},
    {-0.19245008972988, 0, 0.19245008972988, 0.19245008972988, 0,
     -0.19245008972988, -0.19245008972988, 0, 0.19245008972988},
    {0.17023209847088, -0.11111111111111, -0.20882058239687, 0.038588483925985,
     0.22222222222222, 0.038588483925985, -0.20882058239687, -0.11111111111111,
     0.17023209847088},
    {-0.14284169104145, 0.19245008972988, 0.076004476294593, -0.21884616733605,
     0, 0.21884616733605, -0.076004476294593, -0.19245008972988,
     0.14284169104145},
    {0.11111111111111, -0.22222222222222, 0.11111111111111, 0.11111111111111,
     -0.22222222222222, 0.11111111111111, 0.11111111111111, -0.22222222222222,
     0.11111111111111},
    {-0.076004476294593, 0.19245008972988, -0.21884616733605, 0.14284169104145,
     0, -0.14284169104145, 0.21884616733605, -0.19245008972988,
     0.076004476294593},
    {0.038588483925985, -0.11111111111111, 0.17023209847088, -0.20882058239687,
     0.22222222222222, -0.20882058239687, 0.17023209847088, -0.11111111111111,
     0.038588483925985}};

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
  if (i == 0) return 1;
  if (i == 1) return x;
  Real Tm2 = 1;
  Real Tm1 = x;
  Real Ti = 2 * x * Tm1 - Tm2;
  for (int ii = 2; ii <= i; ++ii) {
    Real temp = Ti;
    Ti = 2 * x * Tm1 - Tm2;
    Tm2 = Tm1;
    Tm1 = temp;
  }
  return Ti;
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
