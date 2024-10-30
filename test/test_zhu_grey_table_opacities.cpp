// ======================================================================
// © 2021-2024. Triad National Security, LLC. All rights reserved.  This
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

#include <singularity-opac/photons/opac_photons.hpp>

using namespace singularity;

using pc = PhysicalConstantsCGS;
using DataBox = Spiner::DataBox<Real>;

#ifdef PORTABILITY_STRATEGY_KOKKOS
using atomic_view = Kokkos::MemoryTraits<Kokkos::Atomic>;
#endif

template <typename T>
PORTABLE_INLINE_FUNCTION T FractionalDifference(const T &a, const T &b) {
  return 2 * std::abs(b - a) / (std::abs(a) + std::abs(b) + 1e-20);
}
constexpr Real EPS_TEST = 1e-6;

TEST_CASE("Zhu table photon opacities", "[GrayPhotons]") {
  WHEN("We initialize a gray Zhu tabular photon opacity") {

    // values are directly copied from parsed table for a_grain=-1
    const std::string fbase = "opacitysolar09dustq3p5amax0p1new";
    const std::string fname = fbase + ".hdf5";
    constexpr Real rho_min = 1e-14;  // g/cc.
    constexpr Real temp_min = 1.0; // Kelvin.
    constexpr Real ross_at_min = 0.0030420184427588414; // cm^2/g
    constexpr Real rho_max = 0.7943282347241912;  // g/cc.
    constexpr Real temp_max = 7943282.347242886; // Kelvin.
    constexpr Real ross_at_max = 0.7075511374620657; // cm^2/g

    constexpr Real nu = 3e9; // Hz. UHF microwave

    // not all opacities have a save method, so not using variant for host here
    photons::ZhuTable opac_host = photons::ZhuTable(fname);
    photons::Opacity opac = opac_host.GetOnDevice();

    // Check constants from mean opacity
    THEN("Check constants from mean opacity for consistency") {

      // set test value
      Real mross = opac.AbsorptionCoefficient(rho_min, temp_min, nu);

      // check min rho-T point
      int n_wrong = 0;
      if (FractionalDifference(rho_min * ross_at_min, mross) > EPS_TEST) {
        n_wrong += 1;
      }

      // reset value and check max rho-T point
      mross = opac.AbsorptionCoefficient(rho_max, temp_max, nu);
      if (FractionalDifference(rho_max * ross_at_max, mross) > EPS_TEST) {
        n_wrong += 1;
      }

      REQUIRE(n_wrong == 0);
    }

    // uncomment this save method to creat hdf5 file from ZhuTable DataBox
    // const std::string hdfname = fbase + ".hdf5";
    // opac_host.Save(hdfname);

    // finished
    opac.Finalize();
  }
}
