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

#include <spiner/databox.hpp>
#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>

#include <singularity-opac/base/indexers.hpp>
#include <singularity-opac/base/radiation_types.hpp>
#include <singularity-opac/chebyshev/chebyshev.hpp>
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
constexpr Real EPS_TEST = 1e-3;

TEST_CASE("Gray neutrino opacities", "[GrayNeutrinos]") {
  WHEN("We initialize a gray neutrino opacity") {
    constexpr Real MeV2K = 1e6 * pc::eV / pc::kb;
    constexpr Real MeV2Hz = 1e6 * pc::eV / pc::h;
    constexpr Real rho = 1e11;         // g/cc
    constexpr Real temp = 10 * MeV2K; // 10 MeV
    constexpr Real Ye = 0.1;
    constexpr RadiationType type = RadiationType::NU_ELECTRON;
    constexpr Real nu = 1.25 * MeV2Hz; // 1 MeV

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
            Real jnu = opac.EmissivityPerNuOmega(type, rho, temp, Ye, nu);
            Real Jnu = opac.EmissivityPerNu(type, rho, temp, Ye, nu);
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
      constexpr int nbins = 9;
      constexpr int ntemps = 100;

      Real lnu_min = 0 + std::log10(MeV2Hz);
      Real lnu_max = 1 + std::log10(MeV2Hz);
      Real nu_min = std::pow(10, lnu_min);
      Real nu_max = std::pow(10, lnu_max);

      constexpr Real lt_min = -1;
      constexpr Real lt_max = 2;
      Real dt = (lt_max - lt_min) / (Real)(ntemps - 1);

      Real *nu_bins = (Real *)PORTABLE_MALLOC(nbins * sizeof(Real));
      Real *lnu_bins = (Real *)PORTABLE_MALLOC(nbins * sizeof(Real));
      Real *temp_bins = (Real *)PORTABLE_MALLOC(ntemps * sizeof(Real));
      portableFor(
          "set nu bins", 0, 1, PORTABLE_LAMBDA(const int &i) {
            chebyshev::GetPoints(lnu_min, lnu_max, nbins, lnu_bins);
            for (int j = 0; j < nbins; ++j) {
              nu_bins[j] = std::pow(10, lnu_bins[j]);
            }
          });
      portableFor(
          "set temp bins", 0, ntemps, PORTABLE_LAMBDA(const int &i) {
            temp_bins[i] = std::pow(10, lt_min + dt * i) * MeV2K;
          });

      Real *nu_data = (Real *)PORTABLE_MALLOC(nbins * sizeof(Real));
      Real *lnu_data = (Real *)PORTABLE_MALLOC(nbins * sizeof(Real));
      Real *nu_coeffs = (Real *)PORTABLE_MALLOC(nbins * sizeof(Real));
      
      Real *vm9 = (Real*) PORTABLE_MALLOC(9 * 9 * sizeof(Real));
      portableFor(
          "Fill vm", 0, 1, PORTABLE_LAMBDA(const int& i){ chebyshev::get_vmbox(vm9); });

      portableFor(
          "Fill the indexers", 0, ntemps, PORTABLE_LAMBDA(const int &i) {
            Real temp = temp_bins[i];
            indexers::LogCheb<nbins, Real *> J_cheb(nu_data, lnu_data,
                                                    nu_coeffs, nu_min, nu_max);
            opac.EmissivityPerNu(type, rho, temp, Ye, nu_bins, J_cheb, nbins);
            Real Jtrue = opac.EmissivityPerNu(type, rho, temp, Ye, nu);
            J_cheb.SetInterpCoeffs(Spiner::DataBox(vm9, 9, 9));
            if (std::isnan(J_cheb(nu)) ||
                ((std::abs(Jtrue) >= 1e-14 || J_cheb(nu) >= 1e-14) &&
                 FractionalDifference(J_cheb(nu), Jtrue) > EPS_TEST)) {
              n_wrong_d() += 1;
            }
          });

#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::deep_copy(n_wrong_h, n_wrong_d);
#endif

      REQUIRE(n_wrong_h == 0);

      PORTABLE_FREE(vm9);
      PORTABLE_FREE(nu_data);
      PORTABLE_FREE(lnu_data);
      PORTABLE_FREE(nu_coeffs);
      PORTABLE_FREE(nu_bins);
      PORTABLE_FREE(lnu_bins);
      PORTABLE_FREE(temp_bins);
      
    }

    opac.Finalize();
  }
}

TEST_CASE("Gray photon opacities", "[GrayPhotons]") {
  WHEN("We initialize a gray photon opacity") {
    constexpr Real rho = 1e3;  // g/cc.
    constexpr Real temp = 1e5; // Kelvin.
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
      constexpr int nbins = 100;
      constexpr int ntemps = 100;

      constexpr Real lnu_min = 8;
      constexpr Real lnu_max = 10;
      Real nu_min = std::pow(10, lnu_min);
      Real nu_max = std::pow(10, lnu_max);
      Real dnu = (lnu_max - lnu_min) / (Real)(nbins - 1);

      constexpr Real lt_min = 2;
      constexpr Real lt_max = 4;
      Real dt = (lt_max - lt_min) / (Real)(ntemps - 1);

      Real *nu_bins = (Real *)PORTABLE_MALLOC(nbins * sizeof(Real));
      Real *temp_bins = (Real *)PORTABLE_MALLOC(ntemps * sizeof(Real));
      Real *loglin_bins = (Real *)PORTABLE_MALLOC(ntemps * sizeof(Real));
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
            indexers::LogLinear J_log(loglin_bins, nu_min, nu_max, nbins);
            opac.EmissivityPerNu(rho, temp, nu_bins, J_log, nbins);
            Real Jtrue = opac.EmissivityPerNu(rho, temp, nu);
            if (FractionalDifference(Jtrue, J_log(nu)) > EPS_TEST) {
              n_wrong_d() += 1;
            }
          });

#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::deep_copy(n_wrong_h, n_wrong_d);
#endif
      REQUIRE(n_wrong_h == 0);

      PORTABLE_FREE(nu_bins);
      PORTABLE_FREE(temp_bins);
      PORTABLE_FREE(loglin_bins);
    }

    opac.Finalize();
  }
}

TEST_CASE("Tophat Opacities", "[TopHat]") {
  WHEN("We initialize a tophat neutrino opacity") {
    neutrinos::Opacity opac = neutrinos::Tophat(1, 1e-2, 1e2);
  }
}
