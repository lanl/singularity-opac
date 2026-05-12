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

#ifndef SINGULARITY_OPAC_PHOTONS_POWERLAW_S_OPACITY_PHOTONS_
#define SINGULARITY_OPAC_PHOTONS_POWERLAW_S_OPACITY_PHOTONS_

// This file was partly copied from a file made in part with generative AI.

#include <cassert>
#include <cmath>
#include <cstdio>

#include <ports-of-call/portability.hpp>
#include <singularity-opac/base/opac_error.hpp>

namespace singularity {
namespace photons {

template <typename pc = PhysicalConstantsCGS>
class PowerLawSOpacity {
 public:
  PowerLawSOpacity() = default;
  PowerLawSOpacity(const Real kappa0, const Real rho_exp, const Real temp_exp,
                  const Real nu_exp = 0., const Real nu_ref = 1.,
                  const Real avg_particle_mass = 1.)
      : kappa0_(kappa0), rho_exp_(rho_exp), temp_exp_(temp_exp),
        nu_exp_(nu_exp), nu_ref_(nu_ref), apm_(avg_particle_mass) {
    if (!(nu_ref_ > 0.)) {
      OPAC_ERROR("PowerLawSOpacity: nu_ref must be positive");
    }
  }

  PowerLawSOpacity GetOnDevice() { return *this; }
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return 0; }
  PORTABLE_INLINE_FUNCTION
  void PrintParams() const noexcept {
    printf("Power law scattering opacity. kappa0 = %g rho_exp = %g temp_exp = %g "
           "nu_exp = %g nu_ref = %g avg particle mass = %g\n",
           kappa0_, rho_exp_, temp_exp_, nu_exp_, nu_ref_, apm_);
  }
  inline void Finalize() noexcept {}

  PORTABLE_INLINE_FUNCTION
  Real TotalCrossSection(const Real rho, const Real temp, const Real nu,
                         Real *lambda = nullptr) const {
    return OpacityScale_(rho, temp, nu);
  }

  PORTABLE_INLINE_FUNCTION
  Real DifferentialCrossSection(const Real rho, const Real temp, const Real nu,
                                const Real mu, Real *lambda = nullptr) const {
    // assumed isotropic, elastic
    return OpacityScale_(rho, temp, nu) / (4. * M_PI);
  }

  PORTABLE_INLINE_FUNCTION
  Real TotalScatteringCoefficient(const Real rho, const Real temp,
                                  const Real nu, Real *lambda = nullptr) const {
    return (rho / apm_) * OpacityScale_(rho, temp, nu);
  }

 private:
  PORTABLE_INLINE_FUNCTION
  Real OpacityScale_(const Real rho, const Real temp, const Real nu) const {
    return OpacityPrefactor_(rho, temp) * std::pow(nu / nu_ref_, nu_exp_);
  }

  PORTABLE_INLINE_FUNCTION
  Real OpacityPrefactor_(const Real rho, const Real temp) const {
    return kappa0_ * std::pow(rho, rho_exp_) * std::pow(temp, temp_exp_);
  }

  Real kappa0_;   // Opacity scale. Units depend on nu_exp and nu_ref.
  Real rho_exp_;  // Power law index of density
  Real temp_exp_; // Power law index of temperature
  Real nu_exp_;   // Power law index of frequency
  Real nu_ref_;   // Frequency normalization for nu_exp
  Real apm_;   // Mean molecular weight. Units of g
};

} // namespace photons
} // namespace singularity

#endif // SINGULARITY_OPAC_PHOTONS_POWERLAW_S_OPACITY_PHOTONS_
