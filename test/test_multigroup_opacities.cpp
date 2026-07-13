// ======================================================================
// © 2026. Triad National Security, LLC. All rights reserved.  This
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

// This file was made in part with generative AI.

#include <array>
#include <cmath>
#include <utility>

#include <catch2/catch_test_macros.hpp>

#include <ports-of-call/portability.hpp>
#include <spiner/databox.hpp>

#include <singularity-opac/photons/gray_s_opacity_photons.hpp>
#include <singularity-opac/photons/mean_s_opacity_photons.hpp>
#include <singularity-opac/photons/opac_photons.hpp>

using namespace singularity;

template <typename T>
PORTABLE_INLINE_FUNCTION T FractionalDifference(const T &a, const T &b) {
  return 2 * std::abs(b - a) / (std::abs(a) + std::abs(b) + 1e-100);
}

constexpr Real EPS_TEST = 5e-3;

TEST_CASE("Photon multigroup gray opacities are exact", "[MultigroupPhotons]") {
  constexpr Real rho = 1.;
  constexpr Real temp = 1.e5;
  constexpr Real kappa = 3.e-2;
  constexpr int NRho = 3;
  constexpr int NT = 4;
  constexpr int ngroups = 4;
  constexpr int nnu_per_group = 96;
  const Real lRhoMin = std::log10(0.1 * rho);
  const Real lRhoMax = std::log10(10. * rho);
  const Real lTMin = std::log10(0.25 * temp);
  const Real lTMax = std::log10(4. * temp);
  const std::array<Real, ngroups + 1> group_bounds = {1.e12, 3.e12, 1.e13,
                                                      1.e14, 3.e15};

  photons::Gray opac_host(kappa);
  photons::MeanOpacityBase multigroup_host(
      opac_host, lRhoMin, lRhoMax, NRho, lTMin, lTMax, NT, group_bounds,
      ngroups, nnu_per_group);
  auto multigroup = multigroup_host.GetOnDevice();

  REQUIRE(multigroup_host.ngroups() == ngroups);

  int n_wrong = 0;
  portableReduce(
      "check multigroup gray opacity", 0, ngroups,
      PORTABLE_LAMBDA(const int group, int &accum) {
        const Real alpha_planck =
            multigroup.PlanckGroupAbsorptionCoefficient(rho, temp, group);
        const Real alpha_rosseland =
            multigroup.RosselandGroupAbsorptionCoefficient(rho, temp, group);
        const Real alpha_default =
            multigroup.AbsorptionCoefficient(rho, temp, group);
        if (FractionalDifference(alpha_planck, rho * kappa) > EPS_TEST) {
          accum += 1;
        }
        if (FractionalDifference(alpha_rosseland, rho * kappa) > EPS_TEST) {
          accum += 1;
        }
        if (FractionalDifference(alpha_default, rho * kappa) > EPS_TEST) {
          accum += 1;
        }
      },
      n_wrong);
  REQUIRE(n_wrong == 0);

  multigroup.Finalize();
  multigroup_host.Finalize();
}

TEST_CASE("Photon multigroup can be constructed from pretabulated Spiner data",
          "[MultigroupPhotons]") {
  using DataBox = Spiner::DataBox<Real>;

  constexpr int NRho = 2;
  constexpr int NT = 2;
  constexpr int ngroups = 3;
  constexpr Real kappaP0 = 1.5e-2;
  constexpr Real kappaR0 = 4.0e-3;
  constexpr Real rho_exp_p = 0.25;
  constexpr Real temp_exp_p = -0.5;
  constexpr Real rho_exp_r = -0.2;
  constexpr Real temp_exp_r = 0.3;
  const std::array<Real, ngroups + 1> group_bounds = {2.e11, 3.e11, 4.e11,
                                                      5.e11};
  const Real lRhoMin = -4.;
  const Real lRhoMax = 2.;
  const Real lTMin = 2.;
  const Real lTMax = 8.;

  DataBox kappa_planck(NRho, NT, ngroups);
  kappa_planck.setRange(1, lTMin, lTMax, NT);
  kappa_planck.setRange(2, lRhoMin, lRhoMax, NRho);
  DataBox kappa_rosseland;
  kappa_rosseland.copyMetadata(kappa_planck);

  for (int iRho = 0; iRho < NRho; ++iRho) {
    const Real rho = std::pow(10., kappa_planck.range(2).x(iRho));
    for (int iT = 0; iT < NT; ++iT) {
      const Real temp = std::pow(10., kappa_planck.range(1).x(iT));
      for (int group = 0; group < ngroups; ++group) {
        const Real group_factor = group + 1.;
        kappa_planck(iRho, iT, group) = kappaP0 * std::pow(rho, rho_exp_p) *
                                        std::pow(temp, temp_exp_p) *
                                        group_factor;
        kappa_rosseland(iRho, iT, group) = kappaR0 * std::pow(rho, rho_exp_r) *
                                           std::pow(temp, temp_exp_r) /
                                           group_factor;
      }
    }
  }

  photons::MeanOpacityBase multigroup_host(kappa_planck, kappa_rosseland,
                                                 group_bounds);

  const Real rho_test = std::pow(10., 0.5 * (lRhoMin + lRhoMax));
  const Real temp_test = std::pow(10., 0.5 * (lTMin + lTMax));
  for (int group = 0; group < ngroups; ++group) {
    const Real group_factor = group + 1.;
    const Real kappa_planck_expected = kappaP0 * std::pow(rho_test, rho_exp_p) *
                                       std::pow(temp_test, temp_exp_p) *
                                       group_factor;
    const Real kappa_rosseland_expected =
        kappaR0 * std::pow(rho_test, rho_exp_r) *
        std::pow(temp_test, temp_exp_r) / group_factor;

    REQUIRE(
        FractionalDifference(multigroup_host.PlanckGroupAbsorptionCoefficient(
                                 rho_test, temp_test, group),
                             rho_test * kappa_planck_expected) < EPS_TEST);
    REQUIRE(FractionalDifference(
                multigroup_host.RosselandGroupAbsorptionCoefficient(
                    rho_test, temp_test, group),
                rho_test * kappa_rosseland_expected) < EPS_TEST);
  }

  constexpr Real nu_min = 2.e11;
  constexpr Real nu_max = 4.e11;
  const std::array<Real, ngroups + 1> tail_group_bounds = {
      0., nu_min, nu_max, std::numeric_limits<Real>::infinity()};
  const std::array<Real, ngroups> nu_probe = {
      0.5 * nu_min, std::sqrt(nu_min * nu_max), 2. * nu_max};
  photons::MeanOpacityBase with_tail_bounds(kappa_planck, kappa_rosseland,
                                                  tail_group_bounds);

  REQUIRE(with_tail_bounds.HasGroupBounds());
  REQUIRE(with_tail_bounds.GroupOfNu(0.) == 0);
  REQUIRE(with_tail_bounds.GroupOfNu(0.5 * nu_min) == 0);
  REQUIRE(with_tail_bounds.GroupOfNu(nu_min) == 1);
  REQUIRE(with_tail_bounds.GroupOfNu(nu_max) == 2);
  REQUIRE(with_tail_bounds.GroupOfNu(2. * nu_max) == 2);

  for (int group = 0; group < ngroups; ++group) {
    REQUIRE(FractionalDifference(
                with_tail_bounds.PlanckGroupAbsorptionCoefficientFromNu(
                    rho_test, temp_test, nu_probe[group]),
                multigroup_host.PlanckGroupAbsorptionCoefficient(
                    rho_test, temp_test, group)) < EPS_TEST);
    REQUIRE(FractionalDifference(
                with_tail_bounds.RosselandGroupAbsorptionCoefficientFromNu(
                    rho_test, temp_test, nu_probe[group]),
                multigroup_host.RosselandGroupAbsorptionCoefficient(
                    rho_test, temp_test, group)) < EPS_TEST);
  }

  with_tail_bounds.Finalize();
  multigroup_host.Finalize();
  kappa_planck.finalize();
  kappa_rosseland.finalize();
}

#ifdef SPINER_USE_HDF
TEST_CASE("Photon multigroup tables can round-trip through SP5 HDF",
          "[MultigroupPhotons]") {
  constexpr int NRho = 2;
  constexpr int NT = 2;
  constexpr int ngroups = 3;
  const Real lRhoMin = -4.;
  const Real lRhoMax = 2.;
  const Real lTMin = 2.;
  const Real lTMax = 8.;
  constexpr Real nu_min = 2.e11;
  constexpr Real nu_max = 4.e11;
  const std::array<Real, ngroups + 1> group_bounds = {
      0., nu_min, nu_max, std::numeric_limits<Real>::infinity()};
  const std::array<Real, ngroups> nu_probe = {
      0.5 * nu_min, std::sqrt(nu_min * nu_max), 2. * nu_max};
  using DataBox = Spiner::DataBox<Real>;

  DataBox kappa_planck(NRho, NT, ngroups);
  kappa_planck.setRange(1, lTMin, lTMax, NT);
  kappa_planck.setRange(2, lRhoMin, lRhoMax, NRho);
  DataBox kappa_rosseland;
  kappa_rosseland.copyMetadata(kappa_planck);

  for (int iRho = 0; iRho < NRho; ++iRho) {
    for (int iT = 0; iT < NT; ++iT) {
      for (int group = 0; group < ngroups; ++group) {
        kappa_planck(iRho, iT, group) =
            1.e-2 * (1. + iRho + 2. * iT + 3. * group);
        kappa_rosseland(iRho, iT, group) =
            5.e-3 * (2. + 2. * iRho + iT + group);
      }
    }
  }

  photons::MeanOpacityBase saved(kappa_planck, kappa_rosseland,
                                       group_bounds);
  const char *filename = "multigroup-photon-table.sp5";
  saved.Save(filename);
  photons::MeanOpacityBase loaded(filename);

  REQUIRE(loaded.HasGroupBounds());
  REQUIRE(loaded.ngroups() == ngroups);
  REQUIRE(loaded.GroupOfNu(0.) == 0);
  REQUIRE(loaded.GroupOfNu(0.5 * nu_min) == 0);
  REQUIRE(loaded.GroupOfNu(nu_min) == 1);
  REQUIRE(loaded.GroupOfNu(nu_max) == ngroups - 1);
  REQUIRE(loaded.GroupOfNu(2. * nu_max) == ngroups - 1);

  for (int iRho = 0; iRho < NRho; ++iRho) {
    const Real rho_test =
        std::pow(10., lRhoMin + (lRhoMax - lRhoMin) / (NRho - 1) * iRho);
    for (int iT = 0; iT < NT; ++iT) {
      const Real temp_test =
          std::pow(10., lTMin + (lTMax - lTMin) / (NT - 1) * iT);
      for (int group = 0; group < ngroups; ++group) {
        const Real alpha_planck_expected =
            rho_test * kappa_planck(iRho, iT, group);
        const Real alpha_rosseland_expected =
            rho_test * kappa_rosseland(iRho, iT, group);
        REQUIRE(FractionalDifference(loaded.PlanckGroupAbsorptionCoefficient(
                                         rho_test, temp_test, group),
                                     alpha_planck_expected) < EPS_TEST);
        REQUIRE(FractionalDifference(loaded.RosselandGroupAbsorptionCoefficient(
                                         rho_test, temp_test, group),
                                     alpha_rosseland_expected) < EPS_TEST);
        REQUIRE(
            FractionalDifference(loaded.PlanckGroupAbsorptionCoefficientFromNu(
                                     rho_test, temp_test, nu_probe[group]),
                                 alpha_planck_expected) < EPS_TEST);
        REQUIRE(FractionalDifference(
                    loaded.RosselandGroupAbsorptionCoefficientFromNu(
                        rho_test, temp_test, nu_probe[group]),
                    alpha_rosseland_expected) < EPS_TEST);
      }
    }
  }

  loaded.Finalize();
  saved.Finalize();
  kappa_planck.finalize();
  kappa_rosseland.finalize();
}
#endif

TEST_CASE("Photon multigroup frequency lookup uses half-open group bounds",
          "[MultigroupPhotons]") {
  using DataBox = Spiner::DataBox<Real>;

  constexpr int NRho = 2;
  constexpr int NT = 2;
  constexpr int ngroups = 3;
  constexpr Real rho = 3.;
  constexpr Real temp = 5.e4;
  const Real lRhoMin = std::log10(0.25 * rho);
  const Real lRhoMax = std::log10(4. * rho);
  const Real lTMin = std::log10(0.5 * temp);
  const Real lTMax = std::log10(2. * temp);
  const std::array<Real, ngroups + 1> group_bounds = {1.e12, 2.e12, 4.e12,
                                                      8.e12};

  DataBox kappa_planck(NRho, NT, ngroups);
  kappa_planck.setRange(1, lTMin, lTMax, NT);
  kappa_planck.setRange(2, lRhoMin, lRhoMax, NRho);
  DataBox kappa_rosseland;
  kappa_rosseland.copyMetadata(kappa_planck);

  for (int group = 0; group < ngroups; ++group) {
    for (int iT = 0; iT < NT; ++iT) {
      for (int iRho = 0; iRho < NRho; ++iRho) {
        kappa_planck(iRho, iT, group) = 1.e-2 * (group + 1.);
        kappa_rosseland(iRho, iT, group) = 2.e-2 * (group + 1.);
      }
    }
  }

  photons::MeanOpacityBase multigroup_host(kappa_planck, kappa_rosseland,
                                                 group_bounds);
  photons::MeanOpacity multigroup = multigroup_host;

  REQUIRE(multigroup.HasGroupBounds());
  REQUIRE(multigroup.GroupOfNu(group_bounds[0]) == 0);
  REQUIRE(multigroup.GroupOfNu(1.5e12) == 0);
  REQUIRE(multigroup.GroupOfNu(group_bounds[1]) == 1);
  REQUIRE(multigroup.GroupOfNu(group_bounds[2]) == 2);
  REQUIRE(multigroup.GroupOfNu(group_bounds[ngroups]) == ngroups - 1);

  const Real nu_mid = std::sqrt(group_bounds[1] * group_bounds[2]);
  REQUIRE(
      FractionalDifference(
          multigroup.PlanckGroupAbsorptionCoefficientFromNu(rho, temp, nu_mid),
          multigroup.PlanckGroupAbsorptionCoefficient(rho, temp, 1)) <
      EPS_TEST);
  REQUIRE(FractionalDifference(
              multigroup.RosselandGroupAbsorptionCoefficientFromNu(rho, temp,
                                                                   nu_mid),
              multigroup.RosselandGroupAbsorptionCoefficient(rho, temp, 1)) <
          EPS_TEST);
  REQUIRE(FractionalDifference(
              multigroup.AbsorptionCoefficientFromNu(rho, temp, nu_mid),
              multigroup.AbsorptionCoefficient(rho, temp, 1, photons::Rosseland)) < EPS_TEST);

  multigroup.Finalize();
  kappa_planck.finalize();
  kappa_rosseland.finalize();
}

TEST_CASE("Photon multigroup non-CGS wrapper converts units",
          "[MultigroupPhotons]") {
  using DataBox = Spiner::DataBox<Real>;

  constexpr Real rho = 2.;
  constexpr Real temp = 3.e4;
  constexpr int NRho = 2;
  constexpr int NT = 3;
  constexpr int ngroups = 3;
  constexpr Real time_unit = 13.;
  constexpr Real mass_unit = 17.;
  constexpr Real length_unit = 19.;
  constexpr Real temp_unit = 23.;
  constexpr Real rho_unit =
      mass_unit / (length_unit * length_unit * length_unit);
  const std::array<Real, ngroups + 1> group_bounds = {1.e11, 2.e11, 4.e11,
                                                      8.e11};
  const Real lRhoMin = std::log10(0.25 * rho);
  const Real lRhoMax = std::log10(4. * rho);
  const Real lTMin = std::log10(0.5 * temp);
  const Real lTMax = std::log10(2. * temp);

  DataBox kappa_planck(NRho, NT, ngroups);
  kappa_planck.setRange(1, lTMin, lTMax, NT);
  kappa_planck.setRange(2, lRhoMin, lRhoMax, NRho);
  DataBox kappa_rosseland;
  kappa_rosseland.copyMetadata(kappa_planck);

  for (int i = 0; i < kappa_planck.size(); ++i) {
    kappa_planck(i) = 7.e-3;
    kappa_rosseland(i) = 5.e-3;
  }

  photons::MeanOpacityBase reference_host(kappa_planck, kappa_rosseland,
                                                group_bounds);
  photons::MeanOpacityBase multigroup_base(kappa_planck, kappa_rosseland,
                                                 group_bounds);
  auto funny_host =
      photons::MeanNonCGSUnits<photons::MeanOpacityBase>(
          std::move(multigroup_base), time_unit, mass_unit, length_unit,
          temp_unit);

  REQUIRE(funny_host.ngroups() == ngroups);
  REQUIRE(funny_host.HasGroupBounds());

  for (int group = 0; group < ngroups; ++group) {
    const Real alpha_planck_cgs = funny_host.PlanckGroupAbsorptionCoefficient(
                                      rho / rho_unit, temp / temp_unit, group) /
                                  length_unit;
    const Real alpha_rosseland_cgs =
        funny_host.RosselandGroupAbsorptionCoefficient(
            rho / rho_unit, temp / temp_unit, group) /
        length_unit;

    REQUIRE(
        FractionalDifference(alpha_planck_cgs,
                             reference_host.PlanckGroupAbsorptionCoefficient(
                                 rho, temp, group)) < EPS_TEST);
    REQUIRE(
        FractionalDifference(alpha_rosseland_cgs,
                             reference_host.RosselandGroupAbsorptionCoefficient(
                                 rho, temp, group)) < EPS_TEST);
    REQUIRE(FractionalDifference(funny_host.AbsorptionCoefficient(
                                     rho / rho_unit, temp / temp_unit, group) /
                                     length_unit,
                                 reference_host.AbsorptionCoefficient(
                                     rho, temp, group)) < EPS_TEST);

    const Real nu_mid =
        std::sqrt(group_bounds[group] * group_bounds[group + 1]);
    REQUIRE(funny_host.GroupOfNu(nu_mid * time_unit) == group);
    REQUIRE(FractionalDifference(
                funny_host.PlanckGroupAbsorptionCoefficientFromNu(
                    rho / rho_unit, temp / temp_unit, nu_mid * time_unit) /
                    length_unit,
                reference_host.PlanckGroupAbsorptionCoefficientFromNu(
                    rho, temp, nu_mid)) < EPS_TEST);
    REQUIRE(FractionalDifference(
                funny_host.RosselandGroupAbsorptionCoefficientFromNu(
                    rho / rho_unit, temp / temp_unit, nu_mid * time_unit) /
                    length_unit,
                reference_host.RosselandGroupAbsorptionCoefficientFromNu(
                    rho, temp, nu_mid)) < EPS_TEST);
    REQUIRE(FractionalDifference(
                funny_host.AbsorptionCoefficientFromNu(
                    rho / rho_unit, temp / temp_unit, nu_mid * time_unit) /
                    length_unit,
                reference_host.AbsorptionCoefficientFromNu(rho, temp, nu_mid)) <
            EPS_TEST);
  }

  funny_host.Finalize();
  reference_host.Finalize();
  kappa_planck.finalize();
  kappa_rosseland.finalize();
}

TEST_CASE("Photon multigroup non-CGS wrapper works for monochromatic-built "
          "multigroup opacities",
          "[MultigroupPhotons]") {
  constexpr Real rho = 4.e-6;
  constexpr Real temp = 7.e5;
  constexpr Real kappa = 9.e-3;
  constexpr int NRho = 3;
  constexpr int NT = 4;
  constexpr int ngroups = 4;
  constexpr int nnu_per_group = 96;
  constexpr Real time_unit = 11.;
  constexpr Real mass_unit = 13.;
  constexpr Real length_unit = 17.;
  constexpr Real temp_unit = 19.;
  constexpr Real rho_unit =
      mass_unit / (length_unit * length_unit * length_unit);
  const Real lRhoMin = std::log10(0.25 * rho);
  const Real lRhoMax = std::log10(4. * rho);
  const Real lTMin = std::log10(0.5 * temp);
  const Real lTMax = std::log10(2. * temp);
  const std::array<Real, ngroups + 1> group_bounds = {1.e12, 3.e12, 1.e13,
                                                      1.e14, 3.e15};

  photons::Gray opac_host(kappa);
  photons::MeanOpacityBase reference_host(
      opac_host, lRhoMin, lRhoMax, NRho, lTMin, lTMax, NT, group_bounds,
      ngroups, nnu_per_group);
  photons::MeanOpacityBase multigroup_base(
      opac_host, lRhoMin, lRhoMax, NRho, lTMin, lTMax, NT, group_bounds,
      ngroups, nnu_per_group);
  auto funny_host =
      photons::MeanNonCGSUnits<photons::MeanOpacityBase>(
          std::move(multigroup_base), time_unit, mass_unit, length_unit,
          temp_unit);
  auto funny = funny_host.GetOnDevice();

  REQUIRE(funny.ngroups() == ngroups);
  REQUIRE(funny.HasGroupBounds());

  int n_wrong_h = 0;
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::View<int, atomic_view> n_wrong_d("wrong");
#else
  PortableMDArray<int> n_wrong_d(&n_wrong_h, 1);
#endif

  portableFor(
      "compare multigroup non-cgs wrapper built from monochromatic opacity", 0,
      ngroups, PORTABLE_LAMBDA(const int &group) {
        const Real rho_funny = rho / rho_unit;
        const Real temp_funny = temp / temp_unit;
        const Real alpha_planck_funny = funny.PlanckGroupAbsorptionCoefficient(
            rho_funny, temp_funny, group);
        const Real alpha_rosseland_funny =
            funny.RosselandGroupAbsorptionCoefficient(rho_funny, temp_funny,
                                                      group);
        const Real alpha_default_funny =
            funny.AbsorptionCoefficient(rho_funny, temp_funny, group);
        const Real nu_mid =
            std::sqrt(group_bounds[group] * group_bounds[group + 1]);
        const Real alpha_planck_from_nu_funny =
            funny.PlanckGroupAbsorptionCoefficientFromNu(rho_funny, temp_funny,
                                                         nu_mid * time_unit);
        const Real alpha_rosseland_from_nu_funny =
            funny.RosselandGroupAbsorptionCoefficientFromNu(
                rho_funny, temp_funny, nu_mid * time_unit);
        const Real alpha_default_from_nu_funny =
            funny.AbsorptionCoefficientFromNu(rho_funny, temp_funny,
                                              nu_mid * time_unit);
        const Real alpha_planck =
            reference_host.PlanckGroupAbsorptionCoefficient(rho, temp, group);
        const Real alpha_rosseland =
            reference_host.RosselandGroupAbsorptionCoefficient(rho, temp,
                                                               group);

        if (FractionalDifference(alpha_planck_funny / length_unit,
                                 alpha_planck) > EPS_TEST) {
          n_wrong_d() += 1;
        }
        if (FractionalDifference(alpha_rosseland_funny / length_unit,
                                 alpha_rosseland) > EPS_TEST) {
          n_wrong_d() += 1;
        }
        if (FractionalDifference(alpha_default_funny / length_unit,
                                 alpha_rosseland) > EPS_TEST) {
          n_wrong_d() += 1;
        }
        if (funny.GroupOfNu(nu_mid * time_unit) != group) {
          n_wrong_d() += 1;
        }
        if (FractionalDifference(alpha_planck_from_nu_funny / length_unit,
                                 alpha_planck) > EPS_TEST) {
          n_wrong_d() += 1;
        }
        if (FractionalDifference(alpha_rosseland_from_nu_funny / length_unit,
                                 alpha_rosseland) > EPS_TEST) {
          n_wrong_d() += 1;
        }
        if (FractionalDifference(alpha_default_from_nu_funny / length_unit,
                                 alpha_rosseland) > EPS_TEST) {
          n_wrong_d() += 1;
        }
      });

#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::deep_copy(n_wrong_h, n_wrong_d);
#endif
  REQUIRE(n_wrong_h == 0);

  funny_host.Finalize();
  reference_host.Finalize();
}

TEST_CASE("Photon multigroup with extreme bounds [0, infinity] recovers gray "
          "solution",
          "[MultigroupPhotons]") {
  constexpr Real rho = 3.e-5;
  constexpr Real temp = 2.e6;
  constexpr Real kappa = 5.e-2;
  constexpr int NRho = 4;
  constexpr int NT = 5;
  constexpr int ngroups = 1;
  constexpr int nnu_per_group = 128;
  const Real lRhoMin = std::log10(0.1 * rho);
  const Real lRhoMax = std::log10(10. * rho);
  const Real lTMin = std::log10(0.25 * temp);
  const Real lTMax = std::log10(4. * temp);
  const std::array<Real, ngroups + 1> full_spectrum_bounds = {
      0., std::numeric_limits<Real>::infinity()};

  photons::Gray opac_host(kappa);
  photons::MeanOpacityBase multigroup_host(
      opac_host, lRhoMin, lRhoMax, NRho, lTMin, lTMax, NT, full_spectrum_bounds,
      ngroups, nnu_per_group);

  REQUIRE(multigroup_host.ngroups() == ngroups);
  REQUIRE(multigroup_host.HasGroupBounds());

  // Test at multiple densities and temperatures
  for (int iRho = 0; iRho < NRho; ++iRho) {
    const Real test_rho =
        std::pow(10., lRhoMin + (lRhoMax - lRhoMin) / (NRho - 1) * iRho);
    for (int iT = 0; iT < NT; ++iT) {
      const Real test_temp =
          std::pow(10., lTMin + (lTMax - lTMin) / (NT - 1) * iT);

      const Real alpha_planck =
          multigroup_host.PlanckGroupAbsorptionCoefficient(test_rho, test_temp,
                                                           0);
      const Real alpha_rosseland =
          multigroup_host.RosselandGroupAbsorptionCoefficient(test_rho,
                                                              test_temp, 0);
      const Real alpha_gray = test_rho * kappa;

      // Verify no NaNs or infinities
      REQUIRE(std::isfinite(alpha_planck));
      REQUIRE(std::isfinite(alpha_rosseland));

      // For a gray opacity, both Planck and Rosseland means should equal the
      // gray opacity
      REQUIRE(FractionalDifference(alpha_planck, alpha_gray) < EPS_TEST);
      REQUIRE(FractionalDifference(alpha_rosseland, alpha_gray) < EPS_TEST);
    }
  }

  multigroup_host.Finalize();
}

TEST_CASE("Photon multigroup with tail groups [0, nu_mid] and [nu_mid, "
          "infinity] are numerically stable",
          "[MultigroupPhotons]") {
  constexpr Real rho = 7.e-4;
  constexpr Real temp = 5.e5;
  constexpr Real kappa = 1.2e-2;
  constexpr int NRho = 3;
  constexpr int NT = 4;
  constexpr int ngroups = 2;
  constexpr int nnu_per_group = 96;
  constexpr Real nu_mid = 1.e14;
  const Real lRhoMin = std::log10(0.2 * rho);
  const Real lRhoMax = std::log10(5. * rho);
  const Real lTMin = std::log10(0.5 * temp);
  const Real lTMax = std::log10(2. * temp);
  const std::array<Real, ngroups + 1> tail_bounds = {
      0., nu_mid, std::numeric_limits<Real>::infinity()};

  photons::Gray opac_host(kappa);
  photons::MeanOpacityBase multigroup_host(
      opac_host, lRhoMin, lRhoMax, NRho, lTMin, lTMax, NT, tail_bounds, ngroups,
      nnu_per_group);

  REQUIRE(multigroup_host.ngroups() == ngroups);

  // Test both groups at multiple state points
  for (int iRho = 0; iRho < NRho; ++iRho) {
    const Real test_rho =
        std::pow(10., lRhoMin + (lRhoMax - lRhoMin) / (NRho - 1) * iRho);
    for (int iT = 0; iT < NT; ++iT) {
      const Real test_temp =
          std::pow(10., lTMin + (lTMax - lTMin) / (NT - 1) * iT);

      for (int group = 0; group < ngroups; ++group) {
        const Real alpha_planck =
            multigroup_host.PlanckGroupAbsorptionCoefficient(test_rho,
                                                             test_temp, group);
        const Real alpha_rosseland =
            multigroup_host.RosselandGroupAbsorptionCoefficient(
                test_rho, test_temp, group);
        const Real alpha_default =
            multigroup_host.AbsorptionCoefficient(test_rho, test_temp, group);

        // Verify no NaNs or infinities
        REQUIRE(std::isfinite(alpha_planck));
        REQUIRE(std::isfinite(alpha_rosseland));
        REQUIRE(std::isfinite(alpha_default));

        // Verify positivity
        REQUIRE(alpha_planck >= 0.);
        REQUIRE(alpha_rosseland >= 0.);
        REQUIRE(alpha_default >= 0.);

        // For gray opacity, both groups should give the gray value
        const Real alpha_gray = test_rho * kappa;
        REQUIRE(FractionalDifference(alpha_planck, alpha_gray) < EPS_TEST);
        REQUIRE(FractionalDifference(alpha_rosseland, alpha_gray) < EPS_TEST);
      }
    }
  }

  multigroup_host.Finalize();
}

TEST_CASE("Photon multigroup with asymmetric tail groups is numerically stable",
          "[MultigroupPhotons]") {
  constexpr Real rho = 2.e-3;
  constexpr Real temp = 1.e6;
  constexpr Real kappa = 3.5e-2;
  constexpr int NRho = 3;
  constexpr int NT = 3;
  constexpr int ngroups = 4;
  constexpr int nnu_per_group = 80;
  constexpr Real nu_low = 5.e13;
  constexpr Real nu_high = 2.e15;
  const Real lRhoMin = std::log10(0.25 * rho);
  const Real lRhoMax = std::log10(4. * rho);
  const Real lTMin = std::log10(0.5 * temp);
  const Real lTMax = std::log10(2. * temp);
  // Groups: [0, nu_low], [nu_low, nu_mid], [nu_mid, nu_high], [nu_high, inf]
  const std::array<Real, ngroups + 1> asymmetric_bounds = {
      0., nu_low, std::sqrt(nu_low * nu_high), nu_high,
      std::numeric_limits<Real>::infinity()};

  photons::Gray opac_host(kappa);
  photons::MeanOpacityBase multigroup_host(
      opac_host, lRhoMin, lRhoMax, NRho, lTMin, lTMax, NT, asymmetric_bounds,
      ngroups, nnu_per_group);

  REQUIRE(multigroup_host.ngroups() == ngroups);

  const Real test_rho = std::pow(10., 0.5 * (lRhoMin + lRhoMax));
  const Real test_temp = std::pow(10., 0.5 * (lTMin + lTMax));
  const Real alpha_gray = test_rho * kappa;

  for (int group = 0; group < ngroups; ++group) {
    const Real alpha_planck = multigroup_host.PlanckGroupAbsorptionCoefficient(
        test_rho, test_temp, group);
    const Real alpha_rosseland =
        multigroup_host.RosselandGroupAbsorptionCoefficient(test_rho, test_temp,
                                                            group);

    // Verify no NaNs or infinities
    REQUIRE(std::isfinite(alpha_planck));
    REQUIRE(std::isfinite(alpha_rosseland));

    // Verify positivity
    REQUIRE(alpha_planck >= 0.);
    REQUIRE(alpha_rosseland >= 0.);

    // For gray opacity, all groups should give the gray value
    REQUIRE(FractionalDifference(alpha_planck, alpha_gray) < EPS_TEST);
    REQUIRE(FractionalDifference(alpha_rosseland, alpha_gray) < EPS_TEST);
  }

  multigroup_host.Finalize();
}

TEST_CASE("Photon multigroup with single group [0, nuMax) recovers gray "
          "solution",
          "[MultigroupPhotons]") {
  constexpr Real rho = 4.e-4;
  constexpr Real temp = 8.e5;
  constexpr Real kappa = 7.e-3;
  constexpr int NRho = 3;
  constexpr int NT = 4;
  constexpr int ngroups = 1;
  constexpr int nnu_per_group = 96;
  constexpr Real nu_max = 5.e15;
  const Real lRhoMin = std::log10(0.2 * rho);
  const Real lRhoMax = std::log10(5. * rho);
  const Real lTMin = std::log10(0.25 * temp);
  const Real lTMax = std::log10(4. * temp);
  const std::array<Real, ngroups + 1> low_tail_bounds = {0., nu_max};

  photons::Gray opac_host(kappa);
  photons::MeanOpacityBase multigroup_host(
      opac_host, lRhoMin, lRhoMax, NRho, lTMin, lTMax, NT, low_tail_bounds,
      ngroups, nnu_per_group);

  REQUIRE(multigroup_host.ngroups() == ngroups);
  REQUIRE(multigroup_host.HasGroupBounds());

  // Test at multiple state points
  for (int iRho = 0; iRho < NRho; ++iRho) {
    const Real test_rho =
        std::pow(10., lRhoMin + (lRhoMax - lRhoMin) / (NRho - 1) * iRho);
    for (int iT = 0; iT < NT; ++iT) {
      const Real test_temp =
          std::pow(10., lTMin + (lTMax - lTMin) / (NT - 1) * iT);

      const Real alpha_planck =
          multigroup_host.PlanckGroupAbsorptionCoefficient(test_rho, test_temp,
                                                           0);
      const Real alpha_rosseland =
          multigroup_host.RosselandGroupAbsorptionCoefficient(test_rho,
                                                              test_temp, 0);
      const Real alpha_gray = test_rho * kappa;

      // Verify no NaNs or infinities
      REQUIRE(std::isfinite(alpha_planck));
      REQUIRE(std::isfinite(alpha_rosseland));

      // For gray opacity integrated from [0, nuMax), should equal gray value
      // (nuMax is large enough to capture essentially all Planck weight)
      REQUIRE(FractionalDifference(alpha_planck, alpha_gray) < EPS_TEST);
      REQUIRE(FractionalDifference(alpha_rosseland, alpha_gray) < EPS_TEST);
    }
  }

  // Verify GroupOfNu works correctly
  REQUIRE(multigroup_host.GroupOfNu(0.) == 0);
  REQUIRE(multigroup_host.GroupOfNu(0.5 * nu_max) == 0);
  REQUIRE(multigroup_host.GroupOfNu(nu_max) == 0);

  multigroup_host.Finalize();
}

TEST_CASE("Photon multigroup with single group [nuMin, infinity) recovers "
          "gray solution",
          "[MultigroupPhotons]") {
  constexpr Real rho = 6.e-3;
  constexpr Real temp = 3.e5;
  constexpr Real kappa = 2.5e-2;
  constexpr int NRho = 3;
  constexpr int NT = 4;
  constexpr int ngroups = 1;
  constexpr int nnu_per_group = 128;
  constexpr Real nu_min = 1.e12;
  const Real lRhoMin = std::log10(0.25 * rho);
  const Real lRhoMax = std::log10(4. * rho);
  const Real lTMin = std::log10(0.5 * temp);
  const Real lTMax = std::log10(2. * temp);
  const std::array<Real, ngroups + 1> high_tail_bounds = {
      nu_min, std::numeric_limits<Real>::infinity()};

  photons::Gray opac_host(kappa);
  photons::MeanOpacityBase multigroup_host(
      opac_host, lRhoMin, lRhoMax, NRho, lTMin, lTMax, NT, high_tail_bounds,
      ngroups, nnu_per_group);

  REQUIRE(multigroup_host.ngroups() == ngroups);
  REQUIRE(multigroup_host.HasGroupBounds());

  // Test at multiple state points
  for (int iRho = 0; iRho < NRho; ++iRho) {
    const Real test_rho =
        std::pow(10., lRhoMin + (lRhoMax - lRhoMin) / (NRho - 1) * iRho);
    for (int iT = 0; iT < NT; ++iT) {
      const Real test_temp =
          std::pow(10., lTMin + (lTMax - lTMin) / (NT - 1) * iT);

      const Real alpha_planck =
          multigroup_host.PlanckGroupAbsorptionCoefficient(test_rho, test_temp,
                                                           0);
      const Real alpha_rosseland =
          multigroup_host.RosselandGroupAbsorptionCoefficient(test_rho,
                                                              test_temp, 0);
      const Real alpha_gray = test_rho * kappa;

      // Verify no NaNs or infinities
      REQUIRE(std::isfinite(alpha_planck));
      REQUIRE(std::isfinite(alpha_rosseland));

      // For gray opacity integrated from [nuMin, infinity), should equal gray
      // value (nuMin is small enough to capture essentially all Planck weight)
      REQUIRE(FractionalDifference(alpha_planck, alpha_gray) < EPS_TEST);
      REQUIRE(FractionalDifference(alpha_rosseland, alpha_gray) < EPS_TEST);
    }
  }

  // Verify GroupOfNu works correctly
  REQUIRE(multigroup_host.GroupOfNu(nu_min) == 0);
  REQUIRE(multigroup_host.GroupOfNu(2. * nu_min) == 0);
  REQUIRE(multigroup_host.GroupOfNu(1.e20) == 0);
  REQUIRE(multigroup_host.GroupOfNu(std::numeric_limits<Real>::infinity()) ==
          0);

  multigroup_host.Finalize();
}

TEST_CASE("Photon multigroup GroupOfNu handles extreme bounds correctly",
          "[MultigroupPhotons]") {
  constexpr Real rho = 1.e-2;
  constexpr Real temp = 3.e5;
  constexpr Real kappa = 2.e-2;
  constexpr int NRho = 2;
  constexpr int NT = 2;
  constexpr int ngroups = 3;
  constexpr int nnu_per_group = 64;
  constexpr Real nu_low = 1.e13;
  constexpr Real nu_high = 1.e15;
  const Real lRhoMin = std::log10(0.5 * rho);
  const Real lRhoMax = std::log10(2. * rho);
  const Real lTMin = std::log10(0.5 * temp);
  const Real lTMax = std::log10(2. * temp);
  const std::array<Real, ngroups + 1> extreme_bounds = {
      0., nu_low, nu_high, std::numeric_limits<Real>::infinity()};

  photons::Gray opac_host(kappa);
  photons::MeanOpacityBase multigroup_host(
      opac_host, lRhoMin, lRhoMax, NRho, lTMin, lTMax, NT, extreme_bounds,
      ngroups, nnu_per_group);

  REQUIRE(multigroup_host.ngroups() == ngroups);
  REQUIRE(multigroup_host.HasGroupBounds());

  // Test GroupOfNu at various frequencies
  REQUIRE(multigroup_host.GroupOfNu(0.) == 0);
  REQUIRE(multigroup_host.GroupOfNu(1.e-10) == 0);
  REQUIRE(multigroup_host.GroupOfNu(0.5 * nu_low) == 0);
  REQUIRE(multigroup_host.GroupOfNu(nu_low) == 1);
  REQUIRE(multigroup_host.GroupOfNu(std::sqrt(nu_low * nu_high)) == 1);
  REQUIRE(multigroup_host.GroupOfNu(nu_high) == 2);
  REQUIRE(multigroup_host.GroupOfNu(2. * nu_high) == 2);
  REQUIRE(multigroup_host.GroupOfNu(1.e20) == 2);
  REQUIRE(multigroup_host.GroupOfNu(std::numeric_limits<Real>::infinity()) ==
          ngroups - 1);

  // Verify that AbsorptionCoefficientFromNu works at extreme frequencies
  const Real alpha_low =
      multigroup_host.AbsorptionCoefficientFromNu(rho, temp, 0.);
  const Real alpha_very_low =
      multigroup_host.AbsorptionCoefficientFromNu(rho, temp, 1.e-10);
  const Real alpha_high =
      multigroup_host.AbsorptionCoefficientFromNu(rho, temp, 1.e20);
  const Real alpha_inf = multigroup_host.AbsorptionCoefficientFromNu(
      rho, temp, std::numeric_limits<Real>::infinity());

  // All should be finite and positive
  REQUIRE(std::isfinite(alpha_low));
  REQUIRE(std::isfinite(alpha_very_low));
  REQUIRE(std::isfinite(alpha_high));
  REQUIRE(std::isfinite(alpha_inf));
  REQUIRE(alpha_low > 0.);
  REQUIRE(alpha_very_low > 0.);
  REQUIRE(alpha_high > 0.);
  REQUIRE(alpha_inf > 0.);

  // For gray opacity, all should equal the gray value
  const Real alpha_gray = rho * kappa;
  REQUIRE(FractionalDifference(alpha_low, alpha_gray) < EPS_TEST);
  REQUIRE(FractionalDifference(alpha_very_low, alpha_gray) < EPS_TEST);
  REQUIRE(FractionalDifference(alpha_high, alpha_gray) < EPS_TEST);
  REQUIRE(FractionalDifference(alpha_inf, alpha_gray) < EPS_TEST);

  multigroup_host.Finalize();
}

TEST_CASE("Photon multigroup gray scattering opacities are exact",
          "[MultigroupPhotons][MultigroupScattering]") {
  constexpr Real rho = 1.e-3;
  constexpr Real temp = 1.e6;
  constexpr Real sigma = 0.665e-24; // Thomson cross section
  constexpr Real apm = 2.0e-24;     // avg particle mass ~ 2 protons
  constexpr int NRho = 3;
  constexpr int NT = 4;
  constexpr int ngroups = 4;
  constexpr int nnu_per_group = 96;
  const Real lRhoMin = std::log10(0.1 * rho);
  const Real lRhoMax = std::log10(10. * rho);
  const Real lTMin = std::log10(0.25 * temp);
  const Real lTMax = std::log10(4. * temp);
  const std::array<Real, ngroups + 1> group_bounds = {1.e12, 3.e12, 1.e13,
                                                      1.e14, 3.e15};

  photons::GraySOpacity<PhysicalConstantsCGS> s_opac_host(sigma, apm);
  photons::MeanSOpacityBase multigroup_host(
      s_opac_host, lRhoMin, lRhoMax, NRho, lTMin, lTMax, NT, group_bounds,
      ngroups, nnu_per_group);

  REQUIRE(multigroup_host.ngroups() == ngroups);

  int n_wrong = 0;
  for (int group = 0; group < ngroups; ++group) {
    const Real sigma_planck =
        multigroup_host.PlanckGroupScatteringCoefficient(rho, temp, group);
    const Real sigma_rosseland =
        multigroup_host.RosselandGroupScatteringCoefficient(rho, temp, group);
    const Real sigma_default =
        multigroup_host.ScatteringCoefficient(rho, temp, group);
    const Real sigma_expected = (rho / apm) * sigma;
    if (FractionalDifference(sigma_planck, sigma_expected) > EPS_TEST) {
      n_wrong += 1;
    }
    if (FractionalDifference(sigma_rosseland, sigma_expected) > EPS_TEST) {
      n_wrong += 1;
    }
    if (FractionalDifference(sigma_default, sigma_expected) > EPS_TEST) {
      n_wrong += 1;
    }
  }
  REQUIRE(n_wrong == 0);

  multigroup_host.Finalize();
}

TEST_CASE("Photon multigroup scattering with extreme bounds [0, infinity] "
          "works correctly",
          "[MultigroupPhotons][MultigroupScattering]") {
  constexpr Real rho = 2.e-4;
  constexpr Real temp = 5.e5;
  constexpr Real sigma = 0.665e-24;
  constexpr Real apm = 2.0e-24;
  constexpr int NRho = 3;
  constexpr int NT = 3;
  constexpr int ngroups = 1;
  constexpr int nnu_per_group = 128;
  const Real lRhoMin = std::log10(0.2 * rho);
  const Real lRhoMax = std::log10(5. * rho);
  const Real lTMin = std::log10(0.5 * temp);
  const Real lTMax = std::log10(2. * temp);
  const std::array<Real, ngroups + 1> full_spectrum_bounds = {
      0., std::numeric_limits<Real>::infinity()};

  photons::GraySOpacity<PhysicalConstantsCGS> s_opac_host(sigma, apm);
  photons::MeanSOpacityBase multigroup_host(
      s_opac_host, lRhoMin, lRhoMax, NRho, lTMin, lTMax, NT,
      full_spectrum_bounds, ngroups, nnu_per_group);

  REQUIRE(multigroup_host.ngroups() == ngroups);
  REQUIRE(multigroup_host.HasGroupBounds());

  const Real sigma_expected = (rho / apm) * sigma;
  const Real sigma_planck =
      multigroup_host.PlanckGroupScatteringCoefficient(rho, temp, 0);
  const Real sigma_rosseland =
      multigroup_host.RosselandGroupScatteringCoefficient(rho, temp, 0);

  REQUIRE(std::isfinite(sigma_planck));
  REQUIRE(std::isfinite(sigma_rosseland));
  REQUIRE(FractionalDifference(sigma_planck, sigma_expected) < EPS_TEST);
  REQUIRE(FractionalDifference(sigma_rosseland, sigma_expected) < EPS_TEST);

  multigroup_host.Finalize();
}

TEST_CASE("Photon multigroup scattering GroupOfNu handles extreme bounds",
          "[MultigroupPhotons][MultigroupScattering]") {
  constexpr Real rho = 5.e-3;
  constexpr Real temp = 2.e6;
  constexpr Real sigma = 0.665e-24;
  constexpr Real apm = 2.0e-24;
  constexpr int NRho = 2;
  constexpr int NT = 2;
  constexpr int ngroups = 3;
  constexpr int nnu_per_group = 64;
  constexpr Real nu_low = 1.e13;
  constexpr Real nu_high = 1.e15;
  const Real lRhoMin = std::log10(0.5 * rho);
  const Real lRhoMax = std::log10(2. * rho);
  const Real lTMin = std::log10(0.5 * temp);
  const Real lTMax = std::log10(2. * temp);
  const std::array<Real, ngroups + 1> extreme_bounds = {
      0., nu_low, nu_high, std::numeric_limits<Real>::infinity()};

  photons::GraySOpacity<PhysicalConstantsCGS> s_opac_host(sigma, apm);
  photons::MeanSOpacityBase multigroup_host(
      s_opac_host, lRhoMin, lRhoMax, NRho, lTMin, lTMax, NT, extreme_bounds,
      ngroups, nnu_per_group);

  REQUIRE(multigroup_host.ngroups() == ngroups);
  REQUIRE(multigroup_host.HasGroupBounds());

  // Test GroupOfNu at various frequencies
  REQUIRE(multigroup_host.GroupOfNu(0.) == 0);
  REQUIRE(multigroup_host.GroupOfNu(0.5 * nu_low) == 0);
  REQUIRE(multigroup_host.GroupOfNu(nu_low) == 1);
  REQUIRE(multigroup_host.GroupOfNu(std::sqrt(nu_low * nu_high)) == 1);
  REQUIRE(multigroup_host.GroupOfNu(nu_high) == 2);
  REQUIRE(multigroup_host.GroupOfNu(2. * nu_high) == 2);
  REQUIRE(multigroup_host.GroupOfNu(std::numeric_limits<Real>::infinity()) ==
          ngroups - 1);

  // Verify ScatteringCoefficientFromNu works at extremes
  const Real sigma_expected = (rho / apm) * sigma;
  const Real sigma_low =
      multigroup_host.ScatteringCoefficientFromNu(rho, temp, 0.);
  const Real sigma_high =
      multigroup_host.ScatteringCoefficientFromNu(rho, temp, 1.e20);

  REQUIRE(std::isfinite(sigma_low));
  REQUIRE(std::isfinite(sigma_high));
  REQUIRE(FractionalDifference(sigma_low, sigma_expected) < EPS_TEST);
  REQUIRE(FractionalDifference(sigma_high, sigma_expected) < EPS_TEST);

  multigroup_host.Finalize();
}

#ifdef SPINER_USE_HDF
TEST_CASE("Photon multigroup scattering tables can round-trip through SP5 HDF",
          "[MultigroupPhotons][MultigroupScattering]") {
  using DataBox = Spiner::DataBox<Real>;

  constexpr int NRho = 2;
  constexpr int NT = 2;
  constexpr int ngroups = 3;
  constexpr Real sigma0 = 1.5e-24;
  const Real lRhoMin = -4.;
  const Real lRhoMax = 2.;
  const Real lTMin = 2.;
  const Real lTMax = 8.;
  constexpr Real nu_min = 2.e11;
  constexpr Real nu_max = 4.e11;
  const std::array<Real, ngroups + 1> group_bounds = {
      0., nu_min, nu_max, std::numeric_limits<Real>::infinity()};

  DataBox sigma_planck(NRho, NT, ngroups);
  sigma_planck.setRange(1, lTMin, lTMax, NT);
  sigma_planck.setRange(2, lRhoMin, lRhoMax, NRho);
  DataBox sigma_rosseland;
  sigma_rosseland.copyMetadata(sigma_planck);

  for (int iRho = 0; iRho < NRho; ++iRho) {
    const Real rho = std::pow(10., sigma_planck.range(2).x(iRho));
    for (int iT = 0; iT < NT; ++iT) {
      const Real temp = std::pow(10., sigma_planck.range(1).x(iT));
      for (int group = 0; group < ngroups; ++group) {
        const Real group_factor = group + 1.;
        sigma_planck(iRho, iT, group) = sigma0 * rho * group_factor;
        sigma_rosseland(iRho, iT, group) = sigma0 * rho / group_factor;
      }
    }
  }

  photons::MeanSOpacityBase saved(sigma_planck, sigma_rosseland,
                                        group_bounds);
  const char *filename = "multigroup-photon-scattering-table.sp5";
  saved.Save(filename);
  photons::MeanSOpacityBase loaded(filename);

  REQUIRE(loaded.HasGroupBounds());
  REQUIRE(loaded.ngroups() == ngroups);
  REQUIRE(loaded.GroupOfNu(0.) == 0);
  REQUIRE(loaded.GroupOfNu(nu_min) == 1);
  REQUIRE(loaded.GroupOfNu(nu_max) == ngroups - 1);

  for (int iRho = 0; iRho < NRho; ++iRho) {
    const Real rho_test =
        std::pow(10., lRhoMin + (lRhoMax - lRhoMin) / (NRho - 1) * iRho);
    for (int iT = 0; iT < NT; ++iT) {
      const Real temp_test =
          std::pow(10., lTMin + (lTMax - lTMin) / (NT - 1) * iT);
      for (int group = 0; group < ngroups; ++group) {
        const Real sigma_planck_expected =
            rho_test * sigma_planck(iRho, iT, group) / rho_test;
        const Real sigma_rosseland_expected =
            rho_test * sigma_rosseland(iRho, iT, group) / rho_test;
        REQUIRE(FractionalDifference(loaded.PlanckGroupScatteringCoefficient(
                                         rho_test, temp_test, group),
                                     sigma_planck_expected * rho_test) <
                EPS_TEST);
        REQUIRE(FractionalDifference(loaded.RosselandGroupScatteringCoefficient(
                                         rho_test, temp_test, group),
                                     sigma_rosseland_expected * rho_test) <
                EPS_TEST);
      }
    }
  }

  loaded.Finalize();
  saved.Finalize();
  sigma_planck.finalize();
  sigma_rosseland.finalize();
}
#endif
