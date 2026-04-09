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

#include <catch2/catch.hpp>

#include <ports-of-call/portability.hpp>
#include <spiner/databox.hpp>

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
  photons::MultigroupOpacityBase multigroup_host(
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

TEST_CASE("Photon multigroup can generate logarithmic groups with tails",
          "[MultigroupPhotons]") {
  constexpr Real rho = 1.;
  constexpr Real temp = 1.e5;
  constexpr Real kappa = 3.e-2;
  constexpr int NRho = 3;
  constexpr int NT = 4;
  constexpr int nlog_groups = 2;
  constexpr int ngroups = nlog_groups + 2;
  constexpr int nnu_per_group = 96;
  constexpr Real nu_min = 1.e12;
  constexpr Real nu_max = 1.e16;
  const Real lRhoMin = std::log10(0.1 * rho);
  const Real lRhoMax = std::log10(10. * rho);
  const Real lTMin = std::log10(0.25 * temp);
  const Real lTMax = std::log10(4. * temp);
  const std::array<Real, ngroups> nu_probe = {5.e11, 1.e13, 1.e15, 2.e16};

  photons::Gray opac_host(kappa);
  photons::MultigroupOpacityBase multigroup_host(
      opac_host, lRhoMin, lRhoMax, NRho, lTMin, lTMax, NT, nu_min, nu_max,
      nlog_groups, nnu_per_group);

  REQUIRE(multigroup_host.ngroups() == ngroups);
  REQUIRE(multigroup_host.HasGroupBounds());
  REQUIRE(multigroup_host.GroupOfNu(0.) == 0);
  REQUIRE(multigroup_host.GroupOfNu(0.5 * nu_min) == 0);
  REQUIRE(multigroup_host.GroupOfNu(nu_min) == 1);
  REQUIRE(multigroup_host.GroupOfNu(1.e14) == 2);
  REQUIRE(multigroup_host.GroupOfNu(nu_max) == ngroups - 1);
  REQUIRE(multigroup_host.GroupOfNu(2. * nu_max) == ngroups - 1);

  for (int group = 0; group < ngroups; ++group) {
    REQUIRE(
        FractionalDifference(
            multigroup_host.PlanckGroupAbsorptionCoefficient(rho, temp, group),
            rho * kappa) < EPS_TEST);
    REQUIRE(FractionalDifference(
                multigroup_host.RosselandGroupAbsorptionCoefficient(rho, temp,
                                                                    group),
                rho * kappa) < EPS_TEST);
    REQUIRE(FractionalDifference(multigroup_host.AbsorptionCoefficientFromNu(
                                     rho, temp, nu_probe[group]),
                                 rho * kappa) < EPS_TEST);
  }

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

  photons::MultigroupOpacityBase multigroup_host(kappa_planck, kappa_rosseland,
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
  photons::MultigroupOpacityBase with_tail_bounds(kappa_planck, kappa_rosseland,
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

  photons::MultigroupOpacityBase saved(kappa_planck, kappa_rosseland,
                                       group_bounds);
  const char *filename = "multigroup-photon-table.sp5";
  saved.Save(filename);
  photons::MultigroupOpacityBase loaded(filename);

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

  photons::MultigroupOpacityBase multigroup_host(kappa_planck, kappa_rosseland,
                                                 group_bounds);
  photons::MultigroupOpacity multigroup = multigroup_host;

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
              multigroup.AbsorptionCoefficient(rho, temp, 1)) < EPS_TEST);

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

  photons::MultigroupOpacityBase reference_host(kappa_planck, kappa_rosseland,
                                                group_bounds);
  photons::MultigroupOpacityBase multigroup_base(kappa_planck, kappa_rosseland,
                                                 group_bounds);
  auto funny_host =
      photons::MultigroupNonCGSUnits<photons::MultigroupOpacityBase>(
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
  photons::MultigroupOpacityBase reference_host(
      opac_host, lRhoMin, lRhoMax, NRho, lTMin, lTMax, NT, group_bounds,
      ngroups, nnu_per_group);
  photons::MultigroupOpacityBase multigroup_base(
      opac_host, lRhoMin, lRhoMax, NRho, lTMin, lTMax, NT, group_bounds,
      ngroups, nnu_per_group);
  auto funny_host =
      photons::MultigroupNonCGSUnits<photons::MultigroupOpacityBase>(
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
