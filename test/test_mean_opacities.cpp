// ======================================================================
// © 2022-2024. Triad National Security, LLC. All rights reserved.  This
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
#include <singularity-opac/neutrinos/mean_s_opacity_neutrinos.hpp>
#include <singularity-opac/neutrinos/opac_neutrinos.hpp>
#include <singularity-opac/neutrinos/s_opac_neutrinos.hpp>
#include <singularity-opac/photons/mean_opacity_photons.hpp>
#include <singularity-opac/photons/mean_s_opacity_photons.hpp>
#include <singularity-opac/photons/opac_photons.hpp>
#include <singularity-opac/photons/s_opac_photons.hpp>

using namespace singularity;

using pc = PhysicalConstantsCGS;

#ifdef PORTABILITY_STRATEGY_KOKKOS
using atomic_view = Kokkos::MemoryTraits<Kokkos::Atomic>;
#endif

template <typename T>
PORTABLE_INLINE_FUNCTION T FractionalDifference(const T &a, const T &b) {
  return 2 * std::abs(b - a) / (std::abs(a) + std::abs(b) + 1e-100);
}
constexpr Real EPS_TEST = 1e-3;
template <typename T>
PORTABLE_INLINE_FUNCTION bool IsWrong(const T &a, const T &b) {
  // We only care if the data is different and significantly nonzero.
  constexpr Real ZERO =
      std::numeric_limits<Real>::epsilon() / (EPS_TEST * EPS_TEST);
  return ((std::isnan(a) || std::isnan(b)) ||
          ((std::abs(a) > ZERO || std::abs(b) > ZERO) &&
           FractionalDifference(a, b) > EPS_TEST));
}

TEST_CASE("Mean neutrino opacities", "[MeanNeutrinos]") {
  const std::string grayname = "mean_gray.sp5";

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

    neutrinos::MeanOpacity mean_opac_host = neutrinos::MeanOpacityBase(
        opac_host, lRhoMin, lRhoMax, NRho, lTMin, lTMax, NT, YeMin, YeMax, NYe);
    auto mean_opac = mean_opac_host.GetOnDevice();

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
      neutrinos::MeanOpacity mean_opac_host_load(grayname);
      AND_THEN("The reloaded table matches the gray opacities") {

        auto mean_opac_load = mean_opac_host_load.GetOnDevice();

        int n_wrong = 0;
        portableReduce(
            "rebuilt table vs gray", 0, NRho, 0, NT, 0, NYe, 0, NEUTRINO_NTYPES,
            PORTABLE_LAMBDA(const int iRho, const int iT, const int iYe,
                            const int itp, int &accumulate) {
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

TEST_CASE("Mean neutrino scattering opacities", "[MeanNeutrinosS]") {
  const std::string grayname = "mean_gray_s.sp5";

  WHEN("We initialize a mean neutrino scattering opacity") {
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

    constexpr Real sigma = 1.e-20;

    constexpr Real avg_particle_mass = pc::mp;

    neutrinos::GrayS opac_host(sigma, avg_particle_mass);
    neutrinos::SOpacity opac = opac_host.GetOnDevice();

    neutrinos::MeanSOpacityCGS mean_opac_host(
        opac_host, lRhoMin, lRhoMax, NRho, lTMin, lTMax, NT, YeMin, YeMax, NYe);
    auto mean_opac = mean_opac_host.GetOnDevice();

    THEN("The emissivity per nu omega is consistent with the emissity per nu") {
      int n_wrong_h = 0;
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::View<int, atomic_view> n_wrong_d("wrong");
#else
      PortableMDArray<int> n_wrong_d(&n_wrong_h, 1);
#endif

      portableFor(
          "calc mean opacities", 0, 100, PORTABLE_LAMBDA(const int &i) {
            Real alphaPlanck = mean_opac.PlanckMeanTotalScatteringCoefficient(
                rho, temp, Ye, type);
            Real alphaRosseland =
                mean_opac.PlanckMeanTotalScatteringCoefficient(rho, temp, Ye,
                                                               type);
            if (FractionalDifference(sigma * rho / pc::mp, alphaPlanck) >
                EPS_TEST) {
              n_wrong_d() += 1;
            }
            if (FractionalDifference(sigma * rho / pc::mp, alphaRosseland) >
                EPS_TEST) {
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
      neutrinos::MeanSOpacityCGS mean_opac_host_load(grayname);
      AND_THEN("The reloaded table matches the gray opacities") {

        auto mean_opac_load = mean_opac_host_load.GetOnDevice();

        int n_wrong = 0;
        portableReduce(
            "rebuilt table vs gray", 0, NRho, 0, NT, 0, NYe, 0, NEUTRINO_NTYPES,
            PORTABLE_LAMBDA(const int iRho, const int iT, const int iYe,
                            const int itp, int &accumulate) {
              const Real lRho =
                  lRhoMin + (lRhoMax - lRhoMin) / (NRho - 1) * iRho;
              const Real rho = std::pow(10, lRho);
              const Real lT = lTMin + (lTMax - lTMin) / (NT - 1) * iT;
              const Real T = std::pow(10, lT);
              const Real Ye = YeMin + (YeMax - YeMin) / (NYe - 1) * iYe;
              const RadiationType type = Idx2RadType(itp);

              const Real kappaPgray =
                  mean_opac.PlanckMeanTotalScatteringCoefficient(rho, T, Ye,
                                                                 type);
              const Real kappaPload =
                  mean_opac_load.PlanckMeanTotalScatteringCoefficient(rho, T,
                                                                      Ye, type);

              const Real kappaRgray =
                  mean_opac.RosselandMeanTotalScatteringCoefficient(rho, T, Ye,
                                                                    type);
              const Real kappaRload =
                  mean_opac_load.RosselandMeanTotalScatteringCoefficient(
                      rho, T, Ye, type);

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
          neutrinos::MeanNonCGSUnitsS<neutrinos::MeanSOpacity>(
              std::forward<neutrinos::MeanSOpacity>(mean_opac_host), time_unit,
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
            Real alphaPlanck = mean_opac.PlanckMeanTotalScatteringCoefficient(
                rho, temp, Ye, type);
            Real alphaRosseland =
                mean_opac.RosselandMeanTotalScatteringCoefficient(rho, temp, Ye,
                                                                  type);
            Real alphaPlanckFunny =
                funny_units.PlanckMeanTotalScatteringCoefficient(
                    rho / rho_unit, temp / temp_unit, Ye, type);
            Real alphaRosselandFunny =
                funny_units.RosselandMeanTotalScatteringCoefficient(
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

TEST_CASE("Mean photon opacities", "[MeanPhotons]") {
  const std::string grayname = "mean_gray_photons.sp5";

  WHEN("We initialize a mean photon opacity") {
    constexpr Real MeV2K = 1e6 * pc::eV / pc::kb;
    constexpr Real MeV2Hz = 1e6 * pc::eV / pc::h;
    constexpr Real rho = 1e0;   // g/cc
    constexpr Real temp = 1.e5; // K
    constexpr Real nu = 1.e15;  // Hz

    constexpr int nT = 10;
    constexpr Real lRhoMin = std::log10(0.1 * rho);
    constexpr Real lRhoMax = std::log10(10. * rho);
    constexpr int NRho = 2;
    constexpr Real lTMin = std::log10(0.1 * temp);
    constexpr Real lTMax = std::log10(10. * temp);
    constexpr int NT = 10;

    constexpr Real kappa = 1.e-20;

    photons::Gray opac_host(kappa);
    photons::Opacity opac = opac_host.GetOnDevice();

    photons::MeanOpacity mean_opac_host = photons::MeanOpacityBase(
        opac_host, lRhoMin, lRhoMax, NRho, lTMin, lTMax, NT);
    auto mean_opac = mean_opac_host.GetOnDevice();

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
                mean_opac.PlanckMeanAbsorptionCoefficient(rho, temp);
            Real alphaRosseland =
                mean_opac.PlanckMeanAbsorptionCoefficient(rho, temp);
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
      photons::MeanOpacity mean_opac_host_load =
          photons::MeanOpacityBase(grayname);
      AND_THEN("The reloaded table matches the gray opacities") {

        auto mean_opac_load = mean_opac_host_load.GetOnDevice();

        int n_wrong = 0;
        portableReduce(
            "rebuilt table vs gray", 0, NRho, 0, NT, 0, 0,
            PORTABLE_LAMBDA(const int iRho, const int iT, const int igarbage,
                            int &accumulate) {
              const Real lRho =
                  lRhoMin + (lRhoMax - lRhoMin) / (NRho - 1) * iRho;
              const Real rho = std::pow(10, lRho);
              const Real lT = lTMin + (lTMax - lTMin) / (NT - 1) * iT;
              const Real T = std::pow(10, lT);

              const Real kappaPgray =
                  mean_opac.PlanckMeanAbsorptionCoefficient(rho, T);
              const Real kappaPload =
                  mean_opac_load.PlanckMeanAbsorptionCoefficient(rho, T);

              const Real kappaRgray =
                  mean_opac.RosselandMeanAbsorptionCoefficient(rho, T);
              const Real kappaRload =
                  mean_opac_load.RosselandMeanAbsorptionCoefficient(rho, T);

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

      auto funny_units_host = photons::MeanNonCGSUnits<photons::MeanOpacity>(
          std::forward<photons::MeanOpacity>(mean_opac_host), time_unit,
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
                mean_opac.PlanckMeanAbsorptionCoefficient(rho, temp);
            Real alphaRosseland =
                mean_opac.RosselandMeanAbsorptionCoefficient(rho, temp);
            Real alphaPlanckFunny = funny_units.PlanckMeanAbsorptionCoefficient(
                rho / rho_unit, temp / temp_unit);
            Real alphaRosselandFunny =
                funny_units.RosselandMeanAbsorptionCoefficient(
                    rho / rho_unit, temp / temp_unit);
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

TEST_CASE("Mean photon scattering opacities", "[MeanPhotonS]") {
  const std::string grayname = "mean_gray_s.sp5";

  WHEN("We initialize a mean photon scattering opacity") {
    constexpr Real MeV2K = 1e6 * pc::eV / pc::kb;
    constexpr Real MeV2Hz = 1e6 * pc::eV / pc::h;
    constexpr Real rho = 1e0;   // g/cc
    constexpr Real temp = 1.e5; // K
    constexpr Real nu = 1.e15;  // Hz

    constexpr int nT = 10;
    constexpr Real lRhoMin = std::log10(0.1 * rho);
    constexpr Real lRhoMax = std::log10(10. * rho);
    constexpr int NRho = 2;
    constexpr Real lTMin = std::log10(0.1 * temp);
    constexpr Real lTMax = std::log10(10. * temp);
    constexpr int NT = 10;

    constexpr Real sigma = 1.e-20;

    constexpr Real avg_particle_mass = pc::mp / 2.;

    photons::GrayS opac_host(sigma, avg_particle_mass);
    photons::SOpacity opac = opac_host.GetOnDevice();

    photons::MeanSOpacityCGS mean_opac_host(opac_host, lRhoMin, lRhoMax, NRho,
                                            lTMin, lTMax, NT);
    auto mean_opac = mean_opac_host.GetOnDevice();

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
                mean_opac.PlanckMeanTotalScatteringCoefficient(rho, temp);
            Real alphaRosseland =
                mean_opac.PlanckMeanTotalScatteringCoefficient(rho, temp);
            if (FractionalDifference(sigma * rho / avg_particle_mass,
                                     alphaPlanck) > EPS_TEST) {
              n_wrong_d() += 1;
            }
            if (FractionalDifference(sigma * rho / avg_particle_mass,
                                     alphaRosseland) > EPS_TEST) {
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
      photons::MeanSOpacityCGS mean_opac_host_load(grayname);
      AND_THEN("The reloaded table matches the gray opacities") {

        auto mean_opac_load = mean_opac_host_load.GetOnDevice();

        int n_wrong = 0;
        portableReduce(
            "rebuilt table vs gray", 0, NRho, 0, NT, 0, 0,
            PORTABLE_LAMBDA(const int iRho, const int iT, const int igarbage,
                            int &accumulate) {
              const Real lRho =
                  lRhoMin + (lRhoMax - lRhoMin) / (NRho - 1) * iRho;
              const Real rho = std::pow(10, lRho);
              const Real lT = lTMin + (lTMax - lTMin) / (NT - 1) * iT;
              const Real T = std::pow(10, lT);

              const Real kappaPgray =
                  mean_opac.PlanckMeanTotalScatteringCoefficient(rho, T);
              const Real kappaPload =
                  mean_opac_load.PlanckMeanTotalScatteringCoefficient(rho, T);

              const Real kappaRgray =
                  mean_opac.RosselandMeanTotalScatteringCoefficient(rho, T);
              const Real kappaRload =
                  mean_opac_load.RosselandMeanTotalScatteringCoefficient(rho,
                                                                         T);

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

      auto funny_units_host = photons::MeanNonCGSUnitsS<photons::MeanSOpacity>(
          std::forward<photons::MeanSOpacity>(mean_opac_host), time_unit,
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
                mean_opac.PlanckMeanTotalScatteringCoefficient(rho, temp);
            Real alphaRosseland =
                mean_opac.RosselandMeanTotalScatteringCoefficient(rho, temp);
            Real alphaPlanckFunny =
                funny_units.PlanckMeanTotalScatteringCoefficient(
                    rho / rho_unit, temp / temp_unit);
            Real alphaRosselandFunny =
                funny_units.RosselandMeanTotalScatteringCoefficient(
                    rho / rho_unit, temp / temp_unit);
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
