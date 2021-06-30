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

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <spiner/databox.hpp>

#include <singularity-opac/base/radiation_types.hpp>
#include <singularity-opac/constants/constants.hpp>
#include <singularity-opac/neutrinos/opac_neutrinos.hpp>
#include <singularity-opac/neutrinos/spiner_opac_neutrinos.hpp>

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

TEST_CASE("Spiner opacities, filled with gray data",
          "[GrayNeutrinos][SpinerNeutrinos]") {
  constexpr Real MeV2K = 1e6 * pc::eV / pc::kb;
  constexpr Real lRhoMin = 8;
  constexpr Real lRhoMax = 12;
  constexpr int NRho = 128;
  constexpr Real lTMin = -2*MeV2K;
  constexpr Real lTMax = 2*MeV2K;
  constexpr int NT = 128;
  constexpr Real YeMin = 0.1;
  constexpr Real YeMax = 0.5;
  constexpr int NYe = 64;
  constexpr Real leMin = 1e-1;
  constexpr Real leMax = 1e2;
  constexpr int Ne = 64;
}
