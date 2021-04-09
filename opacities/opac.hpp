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

#ifndef OPACITIES_OPAC_
#define OPACITIES_OPAC_

#ifdef SINGULARITY_ENABLE_EXCEPTIONS
#include <stdexcept>
#define OPAC_ERROR(x) (throw std::runtime_error(x))
#else
#define OPAC_ERROR(x) printf("%s", x)
#endif
#define UNDEFINED_ERROR OPAC_ERROR("DEFINE ME\n")

namespace singularity {

enum class RadiationType { NU_ELECTRON, NU_ANTI, NU_HEAVY, PHOTON };

class ThinFlat {
public:
  ThinFlat(const RadiationType &type) : type_(type) { }
private:
  RadiationType type_;
};

} // namespace singularity

#endif // OPACITIES_OPAC_
