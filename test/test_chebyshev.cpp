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

#include <iostream>

#include <catch2/catch.hpp>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <singularity-opac/chebyshev/chebyshev.hpp>
#include <spiner/databox.hpp>
using namespace singularity::chebyshev;

using DataBox = Spiner::DataBox<Real>;
#ifdef PORTABILITY_STRATEGY_KOKKOS
using atomic_view = Kokkos::MemoryTraits<Kokkos::Atomic>;
#endif

PORTABLE_INLINE_FUNCTION
Real Gauss(Real x, Real mu, Real sigma) {
  return std::exp(-(x - mu) * (x - mu) / (2 * sigma * sigma));
}
template <typename T>
PORTABLE_INLINE_FUNCTION T FractionalDifference(const T &a, const T &b) {
  return 2 * std::abs(b - a) / (std::abs(a) + std::abs(b) + 1e-20);
}
constexpr Real EPS_TEST = 1e-5;

TEST_CASE("Chebyshev Polynomials", "[Chebyshev]") {
  WHEN("We create chebyshev values in x") {
    constexpr int N = 9;
    constexpr Real xmin = 1;
    constexpr Real xmax = 3;
    constexpr Real mu = (xmax + xmin) / 2;
    constexpr Real sigma = 1;
    Real *x = (Real *)PORTABLE_MALLOC(sizeof(Real) * N);
    Real *vm9 = (Real *)PORTABLE_MALLOC(9 * 9 * sizeof(Real));

    portableFor(
        "Set vm", 0, 1, PORTABLE_LAMBDA(const int &i) { get_vmbox(vm9); });

    portableFor(
        "Set x", 0, 1,
        PORTABLE_LAMBDA(const int &i) { GetPoints(xmin, xmax, N, x); });
    WHEN("We can create a Gaussian of x") {
      Real *y = (Real *)PORTABLE_MALLOC(sizeof(Real) * N);
      portableFor(
          "Set y", 0, N,
          PORTABLE_LAMBDA(const int &i) { y[i] = Gauss(x[i], mu, sigma); });
      THEN("We can fit a Chebyshev polynomial") {
        Real *ycoeffs = (Real *)PORTABLE_MALLOC(sizeof(Real) * N);
        portableFor(
            "Compute Cheb polynomial", 0, 1, PORTABLE_LAMBDA(const int &i) {
              MatMultiply(DataBox(vm9, 9, 9), y, ycoeffs, N);
            });
        AND_THEN("The chebyshev polynomials fit") {
          int n_wrong_h = 0;
#ifdef PORTABILITY_STRATEGY_KOKKOS
          Kokkos::View<int, atomic_view> n_wrong_d("wrong");
#else
          PortableMDArray<int> n_wrong_d(&n_wrong_h, 1);
#endif
          constexpr int ntest = 100;
          constexpr Real dx = (xmax - xmin) / (ntest - 1);
          Real *ytest = (Real *)PORTABLE_MALLOC(sizeof(Real) * ntest);
          portableFor(
              "Interpolate", 0, ntest, PORTABLE_LAMBDA(const int &i) {
                Real xtest = xmin + dx * i;
                Real ytest = InterpFromCoeffs(xtest, xmin, xmax, ycoeffs, N);
                Real ytrue = Gauss(xtest, mu, sigma);
                if (std::isnan(ytest) ||
                    FractionalDifference(ytest, ytrue) > EPS_TEST) {
                  n_wrong_d() += 1;
                }
              });
#ifdef PORTABILITY_STRATEGY_KOKKOS
          Kokkos::deep_copy(n_wrong_h, n_wrong_d);
#endif
          REQUIRE(n_wrong_h == 0);
          PORTABLE_FREE(ytest);
        }
        PORTABLE_FREE(ycoeffs);
      }
      PORTABLE_FREE(y);
    }
    PORTABLE_FREE(x);
    PORTABLE_FREE(vm9);
  }
}
