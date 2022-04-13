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
#include <singularity-opac/neutrinos/mean_opacity_neutrinos.hpp>
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

TEST_CASE("Mean neutrino opacities", "[MeanNeutrinos]") {
  WHEN("We initialize a mean neutrino opacity") {
    constexpr Real MeV2K = 1e6 * pc::eV / pc::kb;
    constexpr Real MeV2Hz = 1e6 * pc::eV / pc::h;
    constexpr Real rho = 1e11;        // g/cc
    constexpr Real temp = 10 * MeV2K; // 10 MeV
    constexpr Real Ye = 0.1;
    constexpr RadiationType type = RadiationType::NU_ELECTRON;
    constexpr Real nu = 1.25 * MeV2Hz; // 1 MeV

    constexpr int nT = 10;
    constexpr Real lRhoMin = std::log10(0.1 * rho);
    constexpr Real lRhoMax = std::log10(10. * rho);
    constexpr int NRho = 2;
    constexpr Real lTMin = std::log10(0.1 * temp);
    constexpr Real lTMax = std::log10(10. * temp);
    constexpr int NT = 10;
    constexpr Real YeMin = 0.1;
    constexpr Real YeMax = 0.5;
    constexpr int NYe = 10;

    constexpr Real kappa = 1.e-20;

    neutrinos::Gray opac_host(kappa);
    neutrinos::Opacity opac = opac_host.GetOnDevice();

    neutrinos::MeanOpacity mean_opac_host(opac_host, lRhoMin, lRhoMax, NRho,
                                          lTMin, lTMax, NT, YeMin, YeMax, NYe);
    neutrinos::MeanOpacity mean_opac = mean_opac_host.GetOnDevice();

    THEN("The emissivity per nu omega is consistent with the emissity per nu") {
      int n_wrong_h = 0;
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::View<int, atomic_view> n_wrong_d("wrong");
#else
      PortableMDArray<int> n_wrong_d(&n_wrong_h, 1);
#endif

      portableFor(
          "calc mean opacities", 0, 100, PORTABLE_LAMBDA(const int &i) {
            Real alphaPlanck =
                mean_opac.PlanckMeanAbsorptionCoefficient(rho, temp, Ye, type);
            Real alphaRosseland =
                mean_opac.PlanckMeanAbsorptionCoefficient(rho, temp, Ye, type);
            if (FractionalDifference(kappa * rho, alphaPlanck) > EPS_TEST) {
              n_wrong_d() += 1;
            }
            if (FractionalDifference(kappa * rho, alphaRosseland) > EPS_TEST) {
              n_wrong_d() += 1;
            }
          });

#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::deep_copy(n_wrong_h, n_wrong_d);
#endif
      REQUIRE(n_wrong_h == 0);
    }

#ifdef SPINER_USE_HDF
    THEN("We can save to disk and reload") {
      mean_opac.Save(grayname);
      neutrinos::MeanOpacity mean_opac_host_load =
          neutrinos::MeanOpacity(grayname);
      AND_THEN("The reloaded table matches the gray opacities") {

        neutrinos::MeanOpacity mean_opac_load =
            mean_opac_load_host.GetOnDevice();

        int n_wrong = 0;
        portableReduce(
            "rebuilt table vs gray", 0, NRho, 0, NT, 0, NYe, 0, NEUTRINO_NTYPES,
            0, Ne,
            PORTABLE_LAMBDA(const int iRho, const int iT, const int iYe,
                            const int itp, const int ie, int &accumulate) {
              const Real lRho =
                  lRhoMin + (lRhoMax - lRhoMin) / (NRho - 1) * iRho;
              const Real rho = std::pow(10, lRho);
              const Real lT = lTMin + (lTMax - lTMin) / (NT - 1) * iT;
              const Real T = std::pow(10, lT);
              const Real Ye = YeMin + (YeMax - YeMin) / (NYe - 1) * iYe;
              const RadiationType type = Idx2RadType(itp);

              const Real kappaPgray =
                  mean_opac.PlanckMeanAbsorptionCoefficient(rho, T, Ye, type);
              const Real kappaPload =
                  mean_opac_load.PlanckMeanAbsorptionCoefficient(rho, T, Ye,
                                                                 type);

              const Real kappaRgray =
                  mean_opac.RosselandMeanAbsorptionCoefficient(rho, T, Ye,
                                                               type);
              const Real kappaRload =
                  mean_opac_load.RosselandMeanAbsorptionCoefficient(rho, T, Ye,
                                                                    type);

              if (IsWrong(kappaPgray, kappaPload)) {
                accumulate += 1;
              }
              if (IsWrong(kappaRgray, kappaRload)) {
                accumulate += 1;
              }
            },
            n_wrong);
        REQUIRE(n_wrong == 0);
      }
    }
#endif

    THEN("We can create an opacity object with funny units") {
      constexpr Real time_unit = 123;
      constexpr Real mass_unit = 456;
      constexpr Real length_unit = 789;
      constexpr Real temp_unit = 276;
      constexpr Real rho_unit =
          mass_unit / (length_unit * length_unit * length_unit);

      auto funny_units_host =
          neutrinos::MeanNonCGSUnits<neutrinos::MeanOpacity>(
              std::forward<neutrinos::MeanOpacity>(mean_opac_host), time_unit,
              mass_unit, length_unit, temp_unit);

      auto funny_units = funny_units_host.GetOnDevice();

      int n_wrong_h = 0;
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::View<int, atomic_view> n_wrong_d("wrong");
#else
      PortableMDArray<int> n_wrong_d(&n_wrong_h, 1);
#endif

      portableFor(
          "compare different units", 0, 100, PORTABLE_LAMBDA(const int &i) {
            Real alphaPlanck =
                mean_opac.PlanckMeanAbsorptionCoefficient(rho, temp, Ye, type);
            Real alphaRosseland = mean_opac.RosselandMeanAbsorptionCoefficient(
                rho, temp, Ye, type);
            Real alphaPlanckFunny = funny_units.PlanckMeanAbsorptionCoefficient(
                rho / rho_unit, temp / temp_unit, Ye, type);
            Real alphaRosselandFunny =
                funny_units.RosselandMeanAbsorptionCoefficient(
                    rho / rho_unit, temp / temp_unit, Ye, type);
            if (FractionalDifference(alphaPlanck, alphaPlanckFunny /
                                                      length_unit) > EPS_TEST) {
              n_wrong_d() += 1;
            }
            if (FractionalDifference(alphaRosseland,
                                     alphaRosselandFunny / length_unit) >
                EPS_TEST) {
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
