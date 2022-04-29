// ======================================================================
// © 2022. Triad National Security, LLC. All rights reserved.  This
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
#include <singularity-opac/neutrinos/opac_neutrinos.hpp>
#include <singularity-opac/photons/opac_photons.hpp>

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

TEST_CASE("Neutrino thermal distribution", "[NeutrinoThermalDistribution]") {
  WHEN("We initialize a gray neutrino opacity") {
    constexpr Real MeV2K = 1e6 * pc::eV / pc::kb;
    constexpr Real temp = 10 * MeV2K; // 10 MeV
    constexpr RadiationType type = RadiationType::NU_ELECTRON;

    neutrinos::Gray opac_host(1);
    neutrinos::Opacity opac = opac_host.GetOnDevice();

    THEN("The energy density of temperature is consistent with the temperature "
         "of energy density") {
      int n_wrong_h = 0;
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::View<int, atomic_view> n_wrong_d("wrong");
#else
      PortableMDArray<int> n_wrong_d(&n_wrong_h, 1);
#endif
      portableFor(
          "calc temperatures", 0, 100, PORTABLE_LAMBDA(const int &i) {
            Real Tr0 = (1. + 0.01 * i) * temp;
            Real er = opac.EnergyDensityFromTemperature(Tr0, type);
            Real B = opac.ThermalDistributionOfT(Tr0, type);
            Real Tr = opac.TemperatureFromEnergyDensity(er, type);
            if (FractionalDifference(er, B / pc::c) > EPS_TEST) {
              n_wrong_d() += 1;
            }
            if (FractionalDifference(Tr0, Tr) > EPS_TEST) {
              n_wrong_d() += 1;
            }
          });

#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::deep_copy(n_wrong_h, n_wrong_d);
#endif
      REQUIRE(n_wrong_h == 0);
    }

    opac.Finalize();
  }
}

TEST_CASE("Photon thermal distribution", "[PhotonThermalDistribution]") {
  WHEN("We initialize a gray photon opacity") {
    constexpr Real temp = 1.e6;

    photons::Gray opac_host(1);
    photons::Opacity opac = opac_host.GetOnDevice();

    THEN("The energy density of temperature is consistent with the temperature "
         "of energy density") {
      int n_wrong_h = 0;
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::View<int, atomic_view> n_wrong_d("wrong");
#else
      PortableMDArray<int> n_wrong_d(&n_wrong_h, 1);
#endif
      portableFor(
          "calc temperatures", 0, 100, PORTABLE_LAMBDA(const int &i) {
            Real Tr0 = (1. + 0.01 * i) * temp;
            Real er = opac.EnergyDensityFromTemperature(Tr0);
            Real B = opac.ThermalDistributionOfT(Tr0);
            Real Tr = opac.TemperatureFromEnergyDensity(er);
            if (FractionalDifference(er, B / pc::c) > EPS_TEST) {
              n_wrong_d() += 1;
            }
            if (FractionalDifference(Tr0, Tr) > EPS_TEST) {
              n_wrong_d() += 1;
            }
          });

#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::deep_copy(n_wrong_h, n_wrong_d);
#endif
      REQUIRE(n_wrong_h == 0);
    }

    opac.Finalize();
  }
}
