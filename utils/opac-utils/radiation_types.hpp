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

#ifndef UTILS_OPAC_UTILS_RADIATION_TYPES
#define UTILS_OPAC_UTILS_RADIATION_TYPES

#include "opac_error.hpp"

namespace singularity {

enum class RadiationType {
  TRACER = -1,
  NU_ELECTRON = 0,
  NU_ELECTRON_ANTI = 1,
  NU_HEAVY = 2,
  PHOTON = 3
};

int RadType2Idx(RadiationType type) { return static_cast<int>(type); }
RadiationType Idx2RadType(int i) {
  switch (i) {
  case -1:
    return RadiationType::Tracer;
  case 0:
    return RadiationType::Nu_Electron;
  case 1:
    return RadiationType::NU_ELECTRON_ANTI;
  case 2:
    return RadiationType::NU_HEAVY;
  case 3:
    return RadiationType::PHOTON;
  default:
    OPAC_ERROR("Unknown radiation type");
  }
}

} // namespace singularity

#endif // UTILS_OPAC_UTILS_RADIATION_TYPES
