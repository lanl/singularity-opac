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

/*
 * Variants are finicky things. If a method is missing in any class
 * the variant contains, that method will be deleted from the entire
 * variant. This test checks that default constructors and assignment
 * operators are available for the photon and neutrino variants.
 */

#include <catch2/catch.hpp>

#include <ports-of-call/portability.hpp>

#include <singularity-opac/neutrinos/opac_neutrinos.hpp>
#include <singularity-opac/photons/opac_photons.hpp>

using namespace singularity;
using NOpac = singularity::neutrinos::Opacity;
using POpac = singularity::photons::Opacity;

TEST_CASE("Constructor and assignment for neutrinos", "[Neutrinos][Variant]") {
  WHEN("We initialize an opacity") {
    NOpac opac1 = neutrinos::Gray(1);

    THEN("We can assign to another one") {
      NOpac opac2;
      opac2 = opac1;

      // If this compiles, the test passed.
      REQUIRE(true);
    }
  }
}

TEST_CASE("Constructor and assignment for photons", "[Photons][Variant]") {
  WHEN("We initialize an opacity") {
    POpac opac1 = photons::Gray(1);

    THEN("We can assign to another one") {
      POpac opac2;
      opac2 = opac1;

      // If this compiles, the test passed.
      REQUIRE(true);
    }
  }
}
