// ======================================================================
// © 2022-2026. Triad National Security, LLC. All rights reserved.  This
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

// This file was generated in part with the assistance of generative AI.

#include <array>
#include <cmath>
#include <iostream>

#include <catch2/catch_test_macros.hpp>

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
template <typename Function>
Real IntegrateOverLogNu(const Function &f, const Real nu_min, const Real nu_max,
                        const int nnu) {
  const Real lnu_min = std::log(nu_min);
  const Real lnu_max = std::log(nu_max);
  const Real dlnu = (lnu_max - lnu_min) / (nnu - 1);
  Real integral = 0.;
  for (int inu = 0; inu < nnu; ++inu) {
    const Real weight = (inu == 0 || inu == nnu - 1) ? 0.5 : 1.;
    const Real lnu = lnu_min + inu * dlnu;
    const Real nu = std::exp(lnu);
    integral += weight * f(nu) * nu * dlnu;
  }
  return integral;
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

TEST_CASE("Photon thermal distribution integral scales as T^4",
          "[PhotonThermalDistribution]") {
  constexpr int nnu = 2049;
  const std::array<Real, 3> temps = {1.e4, 3.e5, 2.e6};
  photons::PlanckDistribution<pc> dist;

  std::array<Real, 3> integrated_B{};
  for (int i = 0; i < static_cast<int>(temps.size()); ++i) {
    const Real temp = temps[i];
    const Real nu_scale = pc::kb * temp / pc::h;
    // In x = h nu / (k_B T), these bounds are x in [1e-6, 1e2], which is
    // effectively [0, infinity): the low-frequency tail is negligible as
    // x -> 0, and the high-frequency tail is exponentially suppressed.
    const Real nu_min = 1.e-6 * nu_scale;
    const Real nu_max = 1.e2 * nu_scale;

    const Real Bnu_integral = IntegrateOverLogNu(
        [&](const Real nu) { return dist.ThermalDistributionOfTNu(temp, nu); },
        nu_min, nu_max, nnu);
    integrated_B[i] = 4. * M_PI * Bnu_integral;

    // singularity-opac uses the angle-integrated quantity
    //
    //   ThermalDistributionOfT = 4 pi \int B_nu dnu.
    //
    // The radiation energy density E_r = a T^4. Since
    // E_r = (1 / c) * ThermalDistributionOfT, the corresponding analytic
    // value here is
    //
    //   ThermalDistributionOfT = c a T^4.
    //
    // The per-steradian integral is \int B_nu dnu = (c / 4 pi) a T^4.
    const Real B_from_dist = dist.ThermalDistributionOfT(temp);
    const Real B_analytic = pc::ar * pc::c * std::pow(temp, 4);
    const Real B_analytic_per_sr =
        pc::ar * pc::c * std::pow(temp, 4) / (4. * M_PI);
    REQUIRE(FractionalDifference(integrated_B[i], B_from_dist) < 1.e-8);
    REQUIRE(FractionalDifference(integrated_B[i], B_analytic) < 1.e-8);
    REQUIRE(FractionalDifference(Bnu_integral, B_analytic_per_sr) < 1.e-8);
  }

  const Real expected_ratio = std::pow(temps[2] / temps[0], 4);
  REQUIRE(FractionalDifference(integrated_B[2] / integrated_B[0],
                               expected_ratio) < 1.e-8);
}

TEST_CASE("Photon dBnu/dT matches a finite difference derivative",
          "[PhotonThermalDistribution]") {
  constexpr Real temp = 2.e6;
  constexpr Real rel_step = 1.e-6;
  const std::array<Real, 4> x_values = {1.e-3, 1.e-1, 1., 3.};
  photons::PlanckDistribution<pc> dist;

  for (const Real x : x_values) {
    const Real nu = x * pc::kb * temp / pc::h;
    const Real temp_lo = temp * (1. - rel_step);
    const Real temp_hi = temp * (1. + rel_step);
    const Real dBnudT_fd = (dist.ThermalDistributionOfTNu(temp_hi, nu) -
                            dist.ThermalDistributionOfTNu(temp_lo, nu)) /
                           (temp_hi - temp_lo);
    const Real dBnudT = dist.DThermalDistributionOfTNuDT(temp, nu);

    REQUIRE(FractionalDifference(dBnudT, dBnudT_fd) < 1.e-6);
  }
}
