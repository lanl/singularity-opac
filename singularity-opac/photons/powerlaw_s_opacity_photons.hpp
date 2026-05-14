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
                   const Real nu_exp = 0., const Real nu_ref = 1., const Real nu_off = 0.,
                   const Real rho_ref = 1., const Real rho_off = 0.,
                   const Real temp_ref = 1., const Real temp_off = 0.,
                   const Real avg_particle_mass = 1.)
      : kappa0_(kappa0), rho_exp_(rho_exp), temp_exp_(temp_exp),
        rho_ref_(rho_ref), rho_off_(rho_off), temp_ref_(temp_ref), temp_off_(temp_off),
        nu_exp_(nu_exp), nu_ref_(nu_ref), nu_off_(nu_off), apm_(avg_particle_mass) {
    if (!(nu_ref_ > 0.)) {
      OPAC_ERROR("PowerLawSOpacity: nu_ref must be positive");
    }
    if (!(temp_ref_ > 0.)) {
      OPAC_ERROR("PowerLawSOpacity: temp_ref must be positive");
    }
    if (!(nu_ref_ > 0.)) {
      OPAC_ERROR("PowerLawSOpacity: nu_ref must be positive");
    }
    if (rho_off_ < 0.) {
      OPAC_ERROR("PowerLawSOpacity: rho_off must be nonnegative");
    }
    if (temp_off_ < 0.) {
      OPAC_ERROR("PowerLawSOpacity: temp_off must be nonnegative");
    }
    if (nu_off_ < 0.) {
      OPAC_ERROR("PowerLawSOpacity: nu_off must be nonnegative");
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
    return OpacityPrefactor_(rho, temp) * std::pow((nu + nu_off_) / nu_ref_, nu_exp_);
  }

  PORTABLE_INLINE_FUNCTION
  Real OpacityPrefactor_(const Real rho, const Real temp) const {
    const Real rhom = (rho + rho_off_) / rho_ref_;
    const Real tempm = (temp + temp_off_) / temp_ref_;
    return kappa0_ * std::pow(rhom, rho_exp_) * std::pow(tempm, temp_exp_);
  }

  Real kappa0_;        // Opacity scale. Units depend on nu_exp and nu_ref.
  Real rho_exp_;       // Power law index of density
  Real temp_exp_;      // Power law index of temperature
  Real rho_ref_;       // Density normalization for rho_exp
  Real rho_off_;       // Density offset (same units as rho_ref)
  Real temp_ref_;      // Temperature normalization for temp_exp
  Real temp_off_;      // Temperature offset (same units as temp_ref)
  Real nu_exp_;        // Power law index of frequency
  Real nu_ref_;        // Frequency normalization for nu_exp
  Real nu_off_;        // Frequency offset (same units as nu_ref)
  Real apm_;           // Mean molecular weight. Units of g
};

} // namespace photons
} // namespace singularity

#endif // SINGULARITY_OPAC_PHOTONS_POWERLAW_S_OPACITY_PHOTONS_
