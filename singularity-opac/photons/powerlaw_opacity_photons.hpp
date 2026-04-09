// ======================================================================
// © 2024-2026. Triad National Security, LLC. All rights reserved.  This
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

#ifndef SINGULARITY_OPAC_PHOTONS_POWERLAW_OPACITY_PHOTONS_
#define SINGULARITY_OPAC_PHOTONS_POWERLAW_OPACITY_PHOTONS_

// This file was made in part with generative AI.

#include <cassert>
#include <cmath>
#include <cstdio>

#include <ports-of-call/portability.hpp>
#include <singularity-opac/base/opac_error.hpp>
#include <singularity-opac/photons/thermal_distributions_photons.hpp>

namespace singularity {
namespace photons {

template <typename pc = PhysicalConstantsCGS>
class PowerLawOpacity {
 public:
  using PC = pc;

  PowerLawOpacity() = default;
  PowerLawOpacity(const Real kappa0, const Real rho_exp, const Real temp_exp,
                  const Real nu_exp = 0., const Real nu_ref = 1.)
      : PowerLawOpacity(PlanckDistribution<pc>{}, kappa0, rho_exp, temp_exp,
                        nu_exp, nu_ref) {}
  PowerLawOpacity(const PlanckDistribution<pc> &dist, const Real kappa0,
                  const Real rho_exp, const Real temp_exp, const Real nu_exp,
                  const Real nu_ref = 1.)
      : dist_(dist), kappa0_(kappa0), rho_exp_(rho_exp), temp_exp_(temp_exp),
        nu_exp_(nu_exp), nu_ref_(nu_ref) {
    if (!(nu_ref_ > 0.)) {
      OPAC_ERROR("PowerLawOpacity: nu_ref must be positive");
    }
  }

  PowerLawOpacity GetOnDevice() { return *this; }
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return 0; }
  PORTABLE_INLINE_FUNCTION
  void PrintParams() const noexcept {
    printf("Power law opacity. kappa0 = %g rho_exp = %g temp_exp = %g "
           "nu_exp = %g nu_ref = %g\n",
           kappa0_, rho_exp_, temp_exp_, nu_exp_, nu_ref_);
  }
  inline void Finalize() noexcept {}

  PORTABLE_INLINE_FUNCTION
  Real AbsorptionCoefficient(const Real rho, const Real temp, const Real nu,
                             Real *lambda = nullptr) const {
    return rho * OpacityScale_(rho, temp, nu);
  }

  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void
  AbsorptionCoefficient(const Real rho, const Real temp,
                        FrequencyIndexer &nu_bins, DataIndexer &coeffs,
                        const int nbins, Real *lambda = nullptr) const {
    for (int i = 0; i < nbins; ++i) {
      coeffs[i] = AbsorptionCoefficient(rho, temp, nu_bins[i], lambda);
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real AngleAveragedAbsorptionCoefficient(const Real rho, const Real temp,
                                          const Real nu,
                                          Real *lambda = nullptr) const {
    return AbsorptionCoefficient(rho, temp, nu, lambda);
  }

  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void AngleAveragedAbsorptionCoefficient(
      const Real rho, const Real temp, FrequencyIndexer &nu_bins,
      DataIndexer &coeffs, const int nbins, Real *lambda = nullptr) const {
    for (int i = 0; i < nbins; ++i) {
      coeffs[i] =
          AngleAveragedAbsorptionCoefficient(rho, temp, nu_bins[i], lambda);
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real EmissivityPerNuOmega(const Real rho, const Real temp, const Real nu,
                            Real *lambda = nullptr) const {
    Real Bnu = dist_.ThermalDistributionOfTNu(temp, nu, lambda);
    return rho * OpacityScale_(rho, temp, nu) * Bnu;
  }

  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void
  EmissivityPerNuOmega(const Real rho, const Real temp,
                       FrequencyIndexer &nu_bins, DataIndexer &coeffs,
                       const int nbins, Real *lambda = nullptr) const {
    for (int i = 0; i < nbins; ++i) {
      coeffs[i] = EmissivityPerNuOmega(rho, temp, nu_bins[i], lambda);
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real EmissivityPerNu(const Real rho, const Real temp, const Real nu,
                       Real *lambda = nullptr) const {
    return 4 * M_PI * EmissivityPerNuOmega(rho, temp, nu, lambda);
  }

  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void
  EmissivityPerNu(const Real rho, const Real temp, FrequencyIndexer &nu_bins,
                  DataIndexer &coeffs, const int nbins,
                  Real *lambda = nullptr) const {
    for (int i = 0; i < nbins; ++i) {
      coeffs[i] = EmissivityPerNu(rho, temp, nu_bins[i], lambda);
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real Emissivity(const Real rho, const Real temp,
                  Real *lambda = nullptr) const {
    // Once the opacity depends on frequency, the total emissivity is no longer
    // the gray factorization alpha(T, rho) * ThermalDistributionOfT(T). This
    // class intentionally supports only the monochromatic frequency-resolved
    // emissivity APIs in that case. Supporting the fully integrated form would
    // require carrying factorial/Gamma-style moments and Riemann zeta
    // functions for the frequency-dependent power law.
    if (nu_exp_ != 0.) {
      OPAC_ERROR("PowerLawOpacity: total emissivity is only supported for "
                 "nu_exp = 0. Use EmissivityPerNuOmega or EmissivityPerNu "
                 "for frequency-dependent power laws.");
    }
    return rho * OpacityPrefactor_(rho, temp) *
           dist_.ThermalDistributionOfT(temp, lambda);
  }

  PORTABLE_INLINE_FUNCTION
  Real NumberEmissivity(const Real rho, const Real temp,
                        Real *lambda = nullptr) const {
    // Same limitation as Emissivity(): for frequency-dependent power laws, the
    // integrated number emissivity is not treated as a gray closed-form API,
    // and full support would likewise require factorial/Gamma-style moments
    // and Riemann zeta functions.
    if (nu_exp_ != 0.) {
      OPAC_ERROR("PowerLawOpacity: total number emissivity is only supported "
                 "for nu_exp = 0. Use EmissivityPerNuOmega or EmissivityPerNu "
                 "for frequency-dependent power laws.");
    }
    return OpacityPrefactor_(rho, temp) *
           dist_.ThermalNumberDistributionOfT(temp, lambda);
  }

  PORTABLE_INLINE_FUNCTION
  Real ThermalDistributionOfTNu(const Real temp, const Real nu,
                                Real *lambda = nullptr) const {
    return dist_.ThermalDistributionOfTNu(temp, nu, lambda);
  }

  PORTABLE_INLINE_FUNCTION
  Real DThermalDistributionOfTNuDT(const Real temp, const Real nu,
                                   Real *lambda = nullptr) const {
    return dist_.DThermalDistributionOfTNuDT(temp, nu, lambda);
  }

  PORTABLE_INLINE_FUNCTION
  Real ThermalDistributionOfT(const Real temp, Real *lambda = nullptr) const {
    return dist_.ThermalDistributionOfT(temp, lambda);
  }

  PORTABLE_INLINE_FUNCTION Real
  ThermalNumberDistributionOfT(const Real temp, Real *lambda = nullptr) const {
    return dist_.ThermalNumberDistributionOfT(temp, lambda);
  }

  PORTABLE_INLINE_FUNCTION
  Real EnergyDensityFromTemperature(const Real temp,
                                    Real *lambda = nullptr) const {
    return dist_.EnergyDensityFromTemperature(temp, lambda);
  }

  PORTABLE_INLINE_FUNCTION
  Real TemperatureFromEnergyDensity(const Real er,
                                    Real *lambda = nullptr) const {
    return dist_.TemperatureFromEnergyDensity(er, lambda);
  }

  PORTABLE_INLINE_FUNCTION
  Real NumberDensityFromTemperature(const Real temp,
                                    Real *lambda = nullptr) const {
    return dist_.NumberDensityFromTemperature(temp, lambda);
  }

  PORTABLE_INLINE_FUNCTION RuntimePhysicalConstants
  GetRuntimePhysicalConstants() const {
    return RuntimePhysicalConstants(PC());
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
  PlanckDistribution<pc> dist_;
};

} // namespace photons
} // namespace singularity

#endif // SINGULARITY_OPAC_PHOTONS_POWERLAW_OPACITY_PHOTONS_
