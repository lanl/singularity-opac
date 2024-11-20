// ======================================================================
// © 2024. Triad National Security, LLC. All rights reserved.  This
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

#include <singularity-opac/photons/opac_photons.hpp>

using namespace singularity;

using DataBox = Spiner::DataBox<Real>;

#ifdef PORTABILITY_STRATEGY_KOKKOS
using atomic_view = Kokkos::MemoryTraits<Kokkos::Atomic>;
#endif

template <typename T>
PORTABLE_INLINE_FUNCTION T FractionalDifference(const T &a, const T &b) {
  return 2 * std::abs(b - a) / (std::abs(a) + std::abs(b) + 1e-20);
}
PORTABLE_INLINE_FUNCTION Real CalcFrequency(const int n, const Real nu_min,
                                            const Real nu_max, const int n_nu) {
  const Real lnu_min = std::log(nu_min);
  const Real lnu_max = std::log(nu_max);
  const Real dlnu = (lnu_max - lnu_min) / static_cast<Real>(n_nu);
  return std::exp(lnu_min + (n + 0.5) * dlnu);
}

constexpr Real EPS_TEST = 1e-3;

constexpr Real rho_exp = 2.;
constexpr Real temp_exp = 2.5;

TEST_CASE("Scale free power law photon opacities",
          "[PowerLawScaleFreePhotonOpacities]") {
  WHEN("We initialize a scale-free power law photon opacity") {
    constexpr Real rho = 1.5e0;
    constexpr Real temp = 3.2e0;
    constexpr Real nu_min = 1.e-1;
    constexpr Real nu_max = 1.e1;
    constexpr int n_nu = 100;
    constexpr Real kappa0 = 1.5;

    photons::PowerLawScaleFree opac_host(kappa0, rho_exp, temp_exp);
    photons::Opacity opac = opac_host.GetOnDevice();

    THEN("The emissivity per nu omega is consistent with the emissity per nu") {
      int n_wrong_h = 0;
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::View<int, atomic_view> n_wrong_d("wrong");
#else
      PortableMDArray<int> n_wrong_d(&n_wrong_h, 1);
#endif

      portableFor(
          "calc emissivities", 0, n_nu, PORTABLE_LAMBDA(const int &i) {
            const Real nu = CalcFrequency(i, nu_min, nu_max, n_nu);
            const Real jnu = opac.EmissivityPerNuOmega(rho, temp, nu);
            const Real Jnu = opac.EmissivityPerNu(rho, temp, nu);
            const Real kappa =
                kappa0 * std::pow(rho, rho_exp) * std::pow(temp, temp_exp);
            if (FractionalDifference(jnu, rho * kappa *
                                              opac.ThermalDistributionOfTNu(
                                                  temp, nu)) > EPS_TEST) {
              n_wrong_d() += 1;
            }
            if (FractionalDifference(Jnu, 4 * M_PI * jnu) > EPS_TEST) {
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

TEST_CASE("CGS power law photon opacities", "[PowerLawCGSPhotonOpacities]") {
  constexpr Real rho = 1.5e0;    // g/cc
  constexpr Real temp = 1.e3;    // K
  constexpr Real nu_min = 1.e10; // Hz
  constexpr Real nu_max = 1.e14; // Hz
  constexpr int n_nu = 100;
  constexpr Real kappa0 = 1.5; // cm^2 / g

  WHEN("We initialize a CGS power law photon opacity") {
    photons::PowerLaw opac_host(kappa0, rho_exp, temp_exp);
    photons::Opacity opac = opac_host.GetOnDevice();

    THEN("The emissivity per nu omega is consistent with the emissity per nu") {
      int n_wrong_h = 0;
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::View<int, atomic_view> n_wrong_d("wrong");
#else
      PortableMDArray<int> n_wrong_d(&n_wrong_h, 1);
#endif

      portableFor(
          "calc emissivities", 0, n_nu, PORTABLE_LAMBDA(const int &i) {
            const Real nu = CalcFrequency(i, nu_min, nu_max, n_nu);
            const Real jnu = opac.EmissivityPerNuOmega(rho, temp, nu);
            const Real Jnu = opac.EmissivityPerNu(rho, temp, nu);
            const Real kappa =
                kappa0 * std::pow(rho, rho_exp) * std::pow(temp, temp_exp);
            if (FractionalDifference(jnu, rho * kappa *
                                              opac.ThermalDistributionOfTNu(
                                                  temp, nu)) > EPS_TEST) {
              n_wrong_d() += 1;
            }
            if (FractionalDifference(Jnu, 4 * M_PI * jnu) > EPS_TEST) {
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

  WHEN("We initialize a CGS power law photon opacity with non-CGS units") {
    constexpr Real time_unit = 123.;
    constexpr Real mass_unit = 456.;
    constexpr Real length_unit = 789.;
    constexpr Real temp_unit = 276.;
    constexpr Real rho_unit =
        mass_unit / (length_unit * length_unit * length_unit);
    constexpr Real j_unit = mass_unit / (length_unit * time_unit * time_unit);

    photons::NonCGSUnits<photons::PowerLaw> opac_host(
        photons::PowerLaw(kappa0, rho_exp, temp_exp), time_unit, mass_unit,
        length_unit, temp_unit);
    photons::Opacity opac = opac_host.GetOnDevice();
    photons::PowerLaw opac_cgs_host(kappa0, rho_exp, temp_exp);
    photons::Opacity opac_cgs = opac_cgs_host.GetOnDevice();

    THEN("The emissivity per nu omega is consistent with the emissity per nu") {
      int n_wrong_h = 0;
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::View<int, atomic_view> n_wrong_d("wrong");
#else
      PortableMDArray<int> n_wrong_d(&n_wrong_h, 1);
#endif

      portableFor(
          "calc emissivities", 0, n_nu, PORTABLE_LAMBDA(const int &i) {
            const Real nu = CalcFrequency(i, nu_min, nu_max, n_nu);
            const Real jnu = opac.EmissivityPerNuOmega(
                rho / rho_unit, temp / temp_unit, nu * time_unit);
            const Real Jnu = opac.EmissivityPerNu(
                rho / rho_unit, temp / temp_unit, nu * time_unit);
            const Real kappa =
                kappa0 * std::pow(rho, rho_exp) * std::pow(temp, temp_exp);
            if (FractionalDifference(
                    jnu * j_unit,
                    opac_cgs.EmissivityPerNuOmega(rho, temp, nu)) > EPS_TEST) {
              n_wrong_d() += 1;
            }
            if (FractionalDifference(Jnu, 4 * M_PI * jnu) > EPS_TEST) {
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

