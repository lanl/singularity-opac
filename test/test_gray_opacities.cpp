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

#include <catch2/catch.hpp>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>

#include <singularity-opac/base/indexers.hpp>
#include <singularity-opac/base/radiation_types.hpp>
#include <singularity-opac/constants/constants.hpp>
#include <singularity-opac/neutrinos/opac_neutrinos.hpp>
#include <singularity-opac/photons/opac_photons.hpp>

using namespace singularity;

using pc = PhysicalConstants<CGS>;

#ifdef PORTABILITY_STRATEGY_KOKKOS
using atomic_view = Kokkos::MemoryTraits<Kokkos::Atomic>;
#endif

template <typename T>
PORTABLE_INLINE_FUNCTION T FractionalDifference(const T &a, const T &b) {
  return 2 * std::abs(b - a) / (std::abs(a) + std::abs(b) + 1e-20);
}
constexpr Real EPS_TEST = 1e-5;

TEST_CASE("Gray neutrino opacities", "[GrayNeutrinos]") {
  WHEN("We initialize a gray neutrino opacity") {
    constexpr Real rho = 1e8;                         // g/cc
    constexpr Real temp = 10 * 1e6 * pc::eV / pc::kb; // 10 MeV
    constexpr Real Ye = 0.1;
    constexpr RadiationType type = RadiationType::NU_ELECTRON;
    constexpr Real nu = 1 * 1e6 * pc::eV / pc::h; // 1 MeV

    neutrinos::Opacity opac_host = neutrinos::Gray(1);
    neutrinos::Opacity opac = opac_host.GetOnDevice();
    THEN("The emissivity per nu omega is consistent with the emissity per nu") {
      int n_wrong_h = 0;
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::View<int, atomic_view> n_wrong_d("wrong");
#else
      PortableMDArray<int> n_wrong_d(&n_wrong_h, 1);
#endif

      portableFor(
          "calc emissivities", 0, 100, PORTABLE_LAMBDA(const int &i) {
            Real jnu = opac.EmissivityPerNuOmega(rho, temp, Ye, type, nu);
            Real Jnu = opac.EmissivityPerNu(rho, temp, Ye, type, nu);
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

TEST_CASE("Gray photon opacities", "[GrayPhotons]") {
  WHEN("We initialize a gray photon opacity") {
    constexpr Real rho = 1e2;  // g/cc.
    constexpr Real temp = 300; // Kelvin. Room temperature.
    constexpr Real nu = 3e9;   // Hz. UHF microwave

    photons::Opacity opac_host = photons::Gray(1);
    photons::Opacity opac = opac_host.GetOnDevice();
    THEN("The emissivity per nu omega is consistent with the emissity per nu") {
      int n_wrong_h = 0;
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::View<int, atomic_view> n_wrong_d("wrong");
#else
      PortableMDArray<int> n_wrong_d(&n_wrong_h, 1);
#endif
      portableFor(
          "calc emissivities", 0, 100, PORTABLE_LAMBDA(const int &i) {
            Real jnu = opac.EmissivityPerNuOmega(rho, temp, nu);
            Real Jnu = opac.EmissivityPerNu(rho, temp, nu);
            if (FractionalDifference(Jnu, 4 * M_PI * jnu) > EPS_TEST) {
              n_wrong_d() += 1;
            }
          });
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::deep_copy(n_wrong_h, n_wrong_d);
#endif
      REQUIRE(n_wrong_h == 0);
    }
    THEN("We can fill an indexer allocated in a cell") {
      int n_wrong_h = 0;
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::View<int, atomic_view> n_wrong_d("wrong");
#else
      PortableMDArray<int> n_wrong_d(&n_wrong_h, 1);
#endif
      int nbins = 10;
      int ntemps = 100;

      Real lnu_min = 8;
      Real lnu_max = 10;
      Real nu_min = std::pow(10, lnu_min);
      Real nu_max = std::pow(10, lnu_max);
      Real dnu = (lnu_max - lnu_min) / (Real)(nbins - 1);

      Real lt_min = 2;
      Real lt_max = 4;
      Real dt = (lt_max - lt_min) / (Real)(ntemps - 1);

      Real *nu_bins = (Real *)PORTABLE_MALLOC(nbins * sizeof(Real));
      Real *temp_bins = (Real *)PORTABLE_MALLOC(ntemps * sizeof(Real));
      portableFor(
          "set nu bins", 0, nbins, PORTABLE_LAMBDA(const int &i) {
            nu_bins[i] = std::pow(10, lnu_min + dnu * i);
          });
      portableFor(
          "set temp bins", 0, ntemps, PORTABLE_LAMBDA(const int &i) {
            temp_bins[i] = std::pow(10, lt_min + dt * i);
          });

      portableFor(
          "Fill the indexers", 0, ntemps, PORTABLE_LAMBDA(const int &i) {
            Real temp = temp_bins[i];
            indexers::Linear alpha_lin(nu_min, nu_max, nbins);
            indexers::LogLinear alpha_log(nu_min, nu_max, nbins);
            opac.AbsorptionCoefficientPerNuBin(nu_bins, alpha_lin, nbins, rho,
                                               temp);
            opac.AbsorptionCoefficientPerNuBin(nu_bins, alpha_log, nbins, rho,
                                               temp);
            if (FractionalDifference(alpha_lin(nu), alpha_log(nu)) > EPS_TEST) {
              n_wrong_d() += 1;
            }
          });

#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::deep_copy(n_wrong_h, n_wrong_d);
#endif
      REQUIRE(n_wrong_h == 0);

      PORTABLE_FREE(nu_bins);
      PORTABLE_FREE(temp_bins);
    }
    opac.Finalize();
  }
}

TEST_CASE("Tophat Opacities", "[TopHat]") {
  WHEN("We initialize a tophat neutrino opacity") {
    neutrinos::Opacity opac = neutrinos::Tophat(1, 1e-2, 1e2);
  }
}
