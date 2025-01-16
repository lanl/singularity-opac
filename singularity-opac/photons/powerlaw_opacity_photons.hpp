// ======================================================================
// © 2024. Triad National Security, LLC. All rights reserved.  This
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
  PowerLawOpacity(const Real kappa0, const Real rho_exp, const Real temp_exp)
      : kappa0_(kappa0), rho_exp_(rho_exp), temp_exp_(temp_exp) {}
  PowerLawOpacity(const PlanckDistribution<pc> &dist, const Real kappa0,
                  const Real rho_exp, const Real temp_exp)
      : dist_(dist), kappa0_(kappa0), rho_exp_(rho_exp), temp_exp_(temp_exp) {}

  PowerLawOpacity GetOnDevice() { return *this; }
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return 0; }
  PORTABLE_INLINE_FUNCTION
  void PrintParams() const noexcept {
    printf("Power law opacity. kappa0 = %g rho_exp = %g temp_exp = %g\n",
           kappa0_, rho_exp_, temp_exp_);
  }
  inline void Finalize() noexcept {}

  PORTABLE_INLINE_FUNCTION
  Real AbsorptionCoefficient(const Real rho, const Real temp, const Real nu,
                             Real *lambda = nullptr) const {
    return rho * (kappa0_ * std::pow(rho, rho_exp_) * std::pow(temp, temp_exp_));
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
    return rho *
           (kappa0_ * std::pow(rho, rho_exp_) * std::pow(temp, temp_exp_)) *
           Bnu;
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
    Real B = dist_.ThermalDistributionOfT(temp, lambda);
    return rho *
           (kappa0_ * std::pow(rho, rho_exp_) * std::pow(temp, temp_exp_)) * B;
  }

  PORTABLE_INLINE_FUNCTION
  Real NumberEmissivity(const Real rho, const Real temp,
                        Real *lambda = nullptr) const {
    return (kappa0_ * std::pow(rho, rho_exp_) * std::pow(temp, temp_exp_)) *
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
  Real kappa0_;   // Opacity scale. Units of cm^2/g
  Real rho_exp_;  // Power law index of density
  Real temp_exp_; // Power law index of temperature
  PlanckDistribution<pc> dist_;
};

} // namespace photons
} // namespace singularity

#endif // SINGULARITY_OPAC_PHOTONS_POWERLAW_OPACITY_PHOTONS_
