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

#include <cmath>
#include <iostream>

#include <catch2/catch.hpp>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <spiner/databox.hpp>

#include <singularity-opac/base/indexers.hpp>
#include <singularity-opac/base/radiation_types.hpp>
#include <singularity-opac/chebyshev/chebyshev.hpp>
#include <singularity-opac/constants/constants.hpp>
#include <singularity-opac/neutrinos/s_opac_neutrinos.hpp>
//#include <singularity-opac/photons/opac_photons.hpp>

using namespace singularity;

using pc = PhysicalConstantsCGS;

#ifdef PORTABILITY_STRATEGY_KOKKOS
using atomic_view = Kokkos::MemoryTraits<Kokkos::Atomic>;
#endif

template <typename T>
PORTABLE_INLINE_FUNCTION T FractionalDifference(const T &a, const T &b) {
  return 2 * std::abs(b - a) / (std::abs(a) + std::abs(b) + 1e-20);
}
constexpr Real EPS_TEST = 1e-3;

TEST_CASE("Gray neutrino scattering opacities", "[GraySNeutrinos]") {
  WHEN("We initialize a gray neutrino scatteringopacity") {
    constexpr Real MeV2K = 1e6 * pc::eV / pc::kb;
    constexpr Real MeV2Hz = 1e6 * pc::eV / pc::h;
    constexpr Real rho = 1e11;        // g/cc
    constexpr Real temp = 10 * MeV2K; // 10 MeV
    constexpr Real Ye = 0.1;
    constexpr RadiationType type = RadiationType::NU_ELECTRON;
    constexpr Real nu = 1.25 * MeV2Hz; // 1 MeV

    neutrinos::GrayS opac_host(1);
    neutrinos::SOpacity opac = opac_host.GetOnDevice();

    THEN("The emissivity per nu omega is consistent with the emissity per nu") {
      int n_wrong_h = 0;
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::View<int, atomic_view> n_wrong_d("wrong");
#else
      PortableMDArray<int> n_wrong_d(&n_wrong_h, 1);
#endif

      portableFor(
          "calc s opacities", 0, 100, PORTABLE_LAMBDA(const int &i) {
            Real kappa = opac.ScatteringCoefficient(rho, temp, Ye, type, nu);
            Real kappa_avg = opac.AngleAveragedScatteringCoefficient(
                rho, temp, Ye, type, nu);
            if (FractionalDifference(kappa, kappa_avg) > EPS_TEST) {
              n_wrong_d() += 1;
            }
          });

#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::deep_copy(n_wrong_h, n_wrong_d);
#endif
      REQUIRE(n_wrong_h == 0);
    }

    WHEN("We create a gray opacity object in non-cgs units") {
      constexpr Real time_unit = 123;
      constexpr Real mass_unit = 456;
      constexpr Real length_unit = 789;
      constexpr Real temp_unit = 276;
      constexpr Real rho_unit =
          mass_unit / (length_unit * length_unit * length_unit);
      constexpr Real kappa_unit = 1. / length_unit;
      neutrinos::SOpacity funny_units_host =
          neutrinos::NonCGSUnitsS<neutrinos::GrayS>(neutrinos::GrayS(1),
                                                    time_unit, mass_unit,
                                                    length_unit, temp_unit);
      auto funny_units = funny_units_host.GetOnDevice();

      THEN("We can convert meaningfully into and out of funny units") {
        int n_wrong_h = 0;
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::View<int, atomic_view> n_wrong_d("wrong");
#else
        PortableMDArray<int> n_wrong_d(&n_wrong_h, 1);
#endif

        portableFor(
            "opacities in funny units", 0, 100, PORTABLE_LAMBDA(const int &i) {
              Real kappa_funny = funny_units.ScatteringCoefficient(
                  rho / rho_unit, temp / temp_unit, Ye, type, nu * time_unit);
              Real kappa = opac.ScatteringCoefficient(rho, temp, Ye, type, nu);
              if (FractionalDifference(kappa, kappa_funny * kappa_unit) >
                  EPS_TEST) {
                n_wrong_d() += 1;
              }
            });
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::deep_copy(n_wrong_h, n_wrong_d);
#endif
        REQUIRE(n_wrong_h == 0);
      }
    }
  }
}
