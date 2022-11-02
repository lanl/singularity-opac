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
#include <singularity-opac/photons/s_opac_photons.hpp>

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

TEST_CASE("Thomson photon scattering opacities", "[ThomsonSPhotons]") {
  WHEN("We initialize a Thomson photon scattering opacity") {
    constexpr Real MeV2K = 1e6 * pc::eV / pc::kb;
    constexpr Real MeV2Hz = 1e6 * pc::eV / pc::h;
    constexpr Real rho = 1e0;   // g/cc
    constexpr Real temp = 1.e5; // K
    constexpr Real nu = 1.e15;  // Hz

    constexpr Real avg_particle_mass = pc::mp / 2.;

    photons::ThomsonS opac_host(avg_particle_mass);
    photons::SOpacity opac = opac_host.GetOnDevice();

    const Real sigmaT =
        8. * M_PI / 3. * std::pow(pc::alpha * pc::hbar / (pc::me * pc::c), 2);

    THEN("The emissivity per nu omega is consistent with the emissity per nu") {
      int n_wrong_h = 0;
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::View<int, atomic_view> n_wrong_d("wrong");
#else
      PortableMDArray<int> n_wrong_d(&n_wrong_h, 1);
#endif

      portableFor(
          "calc s opacities", 0, 100, PORTABLE_LAMBDA(const int &i) {
            Real kappa = opac.TotalScatteringCoefficient(rho, temp, nu);
            Real sigma = opac.TotalCrossSection(rho, temp, nu);
            if (FractionalDifference(kappa, (rho / avg_particle_mass) / 2 *
                                                sigma) > EPS_TEST) {
              n_wrong_d() += 1;
            }
            if (FractionalDifference(sigma, sigmaT) > EPS_TEST) {
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
      photons::SOpacity funny_units_host =
          photons::NonCGSUnitsS<photons::ThomsonS>(
              photons::ThomsonS(avg_particle_mass), time_unit, mass_unit,
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
              Real kappa_funny = funny_units.TotalScatteringCoefficient(
                  rho / rho_unit, temp / temp_unit, nu * time_unit);
              Real kappa = opac.TotalScatteringCoefficient(rho, temp, nu);
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
