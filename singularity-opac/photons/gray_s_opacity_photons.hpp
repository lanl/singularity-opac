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

#ifndef SINGULARITY_OPAC_PHOTONS_GRAY_S_OPACITY_PHOTONS_
#define SINGULARITY_OPAC_PHOTONS_GRAY_S_OPACITY_PHOTONS_

#include <cassert>
#include <cmath>
#include <cstdio>

#include <ports-of-call/portability.hpp>
#include <singularity-opac/base/opac_error.hpp>

namespace singularity {
namespace photons {

template <typename pc = PhysicalConstantsCGS>
class GraySOpacity {
 public:
  GraySOpacity() = default;
  GraySOpacity(const Real sigma, const Real avg_particle_mass)
      : sigma_(sigma), apm_(avg_particle_mass) {}

  GraySOpacity GetOnDevice() { return *this; }
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return 0; }
  PORTABLE_INLINE_FUNCTION
  void PrintParams() const noexcept {
    printf("Gray scattering opacity. sigma = %g avg particle mass = %g\n",
           sigma_, apm_);
  }
  inline void Finalize() noexcept {}

  PORTABLE_INLINE_FUNCTION
  Real TotalCrossSection(const Real rho, const Real temp,
                         const Real nu,
                         Real *lambda = nullptr) const {
    return sigma_;
  }

  PORTABLE_INLINE_FUNCTION
  Real DifferentialCrossSection(const Real rho, const Real temp,
                                const Real nu,
                                const Real mu, Real *lambda = nullptr) const {
    return sigma_ / (4. * M_PI);
  }

  PORTABLE_INLINE_FUNCTION
  Real TotalScatteringCoefficient(const Real rho, const Real temp,
                                  const Real nu, Real *lambda = nullptr) const {
    return (rho / apm_) * sigma_;
  }

 private:
  Real sigma_; // Scattering cross section. Units of cm^2
  Real apm_;   // Mean molecular weight. Units of g
};

} // namespace photons
} // namespace singularity

#endif // SINGULARITY_OPAC_PHOTONS_GRAY_S_OPACITY_PHOTONS_
