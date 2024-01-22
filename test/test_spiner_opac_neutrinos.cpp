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

#include <string>

#include <catch2/catch.hpp>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>

#include <singularity-opac/base/indexers.hpp>
#include <singularity-opac/base/radiation_types.hpp>
#include <singularity-opac/constants/constants.hpp>
#include <singularity-opac/neutrinos/opac_neutrinos.hpp>
#include <singularity-opac/neutrinos/spiner_opac_neutrinos.hpp>

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
template <typename T>
PORTABLE_INLINE_FUNCTION bool IsWrong(const T &a, const T &b) {
  // We only care if the data is different and significantly nonzero.
  constexpr Real ZERO =
      std::numeric_limits<Real>::epsilon() / (EPS_TEST * EPS_TEST);
  return ((std::isnan(a) || std::isnan(b)) ||
          ((std::abs(a) > ZERO || std::abs(b) > ZERO) &&
           FractionalDifference(a, b) > EPS_TEST));
}

TEST_CASE("Spiner opacities, filled with gray data",
          "[GrayNeutrinos][SpinerNeutrinos]") {
  constexpr Real MeV2K = 1e6 * pc::eV / pc::kb;
  constexpr Real lRhoMin = 8;
  constexpr Real lRhoMax = 12;
  constexpr int NRho = 32;
  constexpr Real lTMin = -2 + std::log10(MeV2K);
  constexpr Real lTMax = 2 + std::log10(MeV2K);
  constexpr int NT = 32;
  constexpr Real YeMin = 0.1;
  constexpr Real YeMax = 0.5;
  constexpr int NYe = 32;
  constexpr Real leMin = -1;
  constexpr Real leMax = 2;
  constexpr int Ne = 32;
  constexpr Real kappa = 1.0;
  const std::string grayname = "gray.sp5";

  WHEN("We initialize a gray neutrino opacity and tabulate it") {
    using Grid_t = Spiner::RegularGrid1D<Real>;

    // Offset the frequencies to hopefully get over the point
    // where the power laws transition
    Grid_t lRhoGrid(lRhoMin, lRhoMax, NRho);
    Grid_t lTGrid(lTMin, lTMax, NT);
    Grid_t YeGrid(YeMin, YeMax, NYe);
    Grid_t leGrid(leMin, leMax, Ne);

    neutrinos::Opacity gray = neutrinos::Gray(kappa);
    neutrinos::SpinerOpac filled(gray, lRhoMin, lRhoMax, NRho, lTMin, lTMax, NT,
                                 YeMin, YeMax, NYe, leMin, leMax, Ne);

    THEN("The table matches the gray opacities") {
      neutrinos::SpinerOpac opac = filled.GetOnDevice();
      int n_wrong = 0;
      portableReduce(
          "table vs gray depends nu", 0, NRho, 0, NT, 0, NYe, 0,
          NEUTRINO_NTYPES, 0, Ne,
          PORTABLE_LAMBDA(const int iRho, const int iT, const int iYe,
                          const int itp, const int ie, int &accumulate) {
            const Real lRho = lRhoGrid.x(iRho);
            const Real rho = std::pow(10, lRho);
            const Real lT = lTGrid.x(iT);
            const Real T = std::pow(10, lT);
            const Real Ye = YeGrid.x(iYe);
            const Real le = leGrid.x(ie);
            const Real e = std::pow(10, le);
            const Real nu = e * neutrinos::SpinerOpac::MeV2Hz;
            const RadiationType type = Idx2RadType(itp);

            // alphanu
            const Real alpha_gray =
                gray.AbsorptionCoefficient(rho, T, Ye, type, nu);
            const Real alpha_tabl =
                opac.AbsorptionCoefficient(rho, T, Ye, type, nu);
            if (IsWrong(alpha_gray, alpha_tabl)) {
              accumulate += 1;
            }

            // Alphanu
            const Real Alpha_gray =
                gray.AngleAveragedAbsorptionCoefficient(rho, T, Ye, type, nu);
            const Real Alpha_tabl =
                opac.AngleAveragedAbsorptionCoefficient(rho, T, Ye, type, nu);
            if (IsWrong(Alpha_gray, Alpha_tabl)) {
              accumulate += 1;
            }

            // jnu
            const Real jgray = gray.EmissivityPerNuOmega(rho, T, Ye, type, nu);
            const Real jtable = opac.EmissivityPerNuOmega(rho, T, Ye, type, nu);
            if (IsWrong(jgray, jtable)) {
              accumulate += 1;
            }
          },
          n_wrong);
      REQUIRE(n_wrong == 0);

      Real *nu_bins = (Real *)PORTABLE_MALLOC(Ne * sizeof(Real));
      portableFor(
          "fill nu bins", 0, Ne, PORTABLE_LAMBDA(const int ie) {
            const Real le = leGrid.x(ie);
            const Real e = std::pow(10, le);
            nu_bins[ie] = e * neutrinos::SpinerOpac::MeV2Hz;
          });

      n_wrong = 0;
      portableReduce(
          "table vs gray indexer API", 0, NRho, 0, NT, 0, NYe, 0,
          NEUTRINO_NTYPES,
          PORTABLE_LAMBDA(const int iRho, const int iT, const int iYe,
                          const int itp, int &accumulate) {
            const Real lRho = lRhoGrid.x(iRho);
            const Real rho = std::pow(10, lRho);
            const Real lT = lTGrid.x(iT);
            const Real T = std::pow(10, lT);
            const Real Ye = YeGrid.x(iYe);
            const RadiationType type = Idx2RadType(itp);
            Real data_gray[Ne];
            Real data_tabl[Ne];

            // alphanu
            gray.AbsorptionCoefficient(rho, T, Ye, type, nu_bins, data_gray,
                                       Ne);
            opac.AbsorptionCoefficient(rho, T, Ye, type, nu_bins, data_tabl,
                                       Ne);
            for (int ie = 0; ie < Ne; ++ie) {
              if (IsWrong(data_gray[ie], data_tabl[ie])) {
                accumulate += 1;
              }
            }

            // Alphanu
            gray.AngleAveragedAbsorptionCoefficient(rho, T, Ye, type, nu_bins,
                                                    data_gray, Ne);
            opac.AngleAveragedAbsorptionCoefficient(rho, T, Ye, type, nu_bins,
                                                    data_tabl, Ne);
            for (int ie = 0; ie < Ne; ++ie) {
              if (IsWrong(data_gray[ie], data_tabl[ie])) {
                accumulate += 1;
              }
            }

            // jnu
            gray.EmissivityPerNuOmega(rho, T, Ye, type, nu_bins, data_gray, Ne);
            opac.EmissivityPerNuOmega(rho, T, Ye, type, nu_bins, data_tabl, Ne);
            for (int ie = 0; ie < Ne; ++ie) {
              if (IsWrong(data_gray[ie], data_tabl[ie])) {
                accumulate += 1;
              }
            }

            // Jnu
            gray.EmissivityPerNu(rho, T, Ye, type, nu_bins, data_gray, Ne);
            opac.EmissivityPerNu(rho, T, Ye, type, nu_bins, data_tabl, Ne);
            for (int ie = 0; ie < Ne; ++ie) {
              if (IsWrong(data_gray[ie], data_tabl[ie])) {
                accumulate += 1;
              }
            }
          },
          n_wrong);
      PORTABLE_FREE(nu_bins);
      REQUIRE(n_wrong == 0);

      n_wrong = 0;
      portableReduce(
          "table vs gray totals", 0, NRho, 0, NT, 0, NYe, 0, NEUTRINO_NTYPES,
          PORTABLE_LAMBDA(const int iRho, const int iT, const int iYe,
                          const int itp, int &accumulate) {
            const Real lRho = lRhoGrid.x(iRho);
            const Real rho = std::pow(10, lRho);
            const Real lT = lTGrid.x(iT);
            const Real T = std::pow(10, lT);
            const Real Ye = YeGrid.x(iYe);
            const RadiationType type = Idx2RadType(itp);

            // Jnu
            const Real Jgray = gray.Emissivity(rho, T, Ye, type);
            const Real Jtabl = opac.Emissivity(rho, T, Ye, type);
            if (IsWrong(Jgray, Jtabl)) {
              accumulate += 1;
            }

            // JYe
            const Real JYe_gray = gray.NumberEmissivity(rho, T, Ye, type);
            const Real JYe_tabl = opac.NumberEmissivity(rho, T, Ye, type);
            if (IsWrong(JYe_gray, JYe_tabl)) {
              accumulate += 1;
            }
          },
          n_wrong);
      REQUIRE(n_wrong == 0);

      opac.Finalize();
    }

#ifdef SPINER_USE_HDF
    THEN("We can save to disk and reload") {
      filled.Save(grayname);
      neutrinos::Opacity opac_host = neutrinos::SpinerOpac(grayname);
      AND_THEN("The reloaded table matches the gray opacities") {

        neutrinos::Opacity opac = opac_host.GetOnDevice();

        int n_wrong = 0;
        portableReduce(
            "rebuilt table vs gray", 0, NRho, 0, NT, 0, NYe, 0, NEUTRINO_NTYPES,
            0, Ne,
            PORTABLE_LAMBDA(const int iRho, const int iT, const int iYe,
                            const int itp, const int ie, int &accumulate) {
              const Real lRho = lRhoGrid.x(iRho);
              const Real rho = std::pow(10, lRho);
              const Real lT = lTGrid.x(iT);
              const Real T = std::pow(10, lT);
              const Real Ye = YeGrid.x(iYe);
              const Real le = leGrid.x(ie);
              const Real e = std::pow(10, le);
              const Real nu = e * neutrinos::SpinerOpac::MeV2Hz;
              const RadiationType type = Idx2RadType(itp);
              const Real Jgray =
                  gray.EmissivityPerNuOmega(rho, T, Ye, type, nu);
              const Real Jtable =
                  filled.EmissivityPerNuOmega(rho, T, Ye, type, nu);

              // alphanu
              const Real alpha_gray =
                  gray.AbsorptionCoefficient(rho, T, Ye, type, nu);
              const Real alpha_tabl =
                  opac.AbsorptionCoefficient(rho, T, Ye, type, nu);
              if (IsWrong(alpha_gray, alpha_tabl)) {
                accumulate += 1;
              }

              // Alpha
              const Real Alpha_gray =
                  gray.AngleAveragedAbsorptionCoefficient(rho, T, Ye, type, nu);
              const Real Alpha_tabl =
                  opac.AngleAveragedAbsorptionCoefficient(rho, T, Ye, type, nu);
              if (IsWrong(Alpha_gray, Alpha_tabl)) {
                accumulate += 1;
              }

              // jnu
              const Real jgray =
                  gray.EmissivityPerNuOmega(rho, T, Ye, type, nu);
              const Real jtable =
                  opac.EmissivityPerNuOmega(rho, T, Ye, type, nu);
              if (IsWrong(jgray, jtable)) {
                accumulate += 1;
              }
            },
            n_wrong);
        REQUIRE(n_wrong == 0);

        n_wrong = 0;
        portableReduce(
            "table vs gray totals", 0, NRho, 0, NT, 0, NYe, 0, NEUTRINO_NTYPES,
            PORTABLE_LAMBDA(const int iRho, const int iT, const int iYe,
                            const int itp, int &accumulate) {
              const Real lRho = lRhoGrid.x(iRho);
              const Real rho = std::pow(10, lRho);
              const Real lT = lTGrid.x(iT);
              const Real T = std::pow(10, lT);
              const Real Ye = YeGrid.x(iYe);
              const RadiationType type = Idx2RadType(itp);

              // Jnu
              const Real Jgray = gray.Emissivity(rho, T, Ye, type);
              const Real Jtabl = opac.Emissivity(rho, T, Ye, type);
              if (IsWrong(Jgray, Jtabl)) {
                accumulate += 1;
              }

              // JYe
              const Real JYe_gray = gray.NumberEmissivity(rho, T, Ye, type);
              const Real JYe_tabl = opac.NumberEmissivity(rho, T, Ye, type);
              if (IsWrong(JYe_gray, JYe_tabl)) {
                accumulate += 1;
              }
            },
            n_wrong);
        REQUIRE(n_wrong == 0);

        // Only need to finalize this one. Gray has nothing to
        // finalize
        opac.Finalize();
      }
    }
#endif // SPINER_USE_HDF
  }
}
