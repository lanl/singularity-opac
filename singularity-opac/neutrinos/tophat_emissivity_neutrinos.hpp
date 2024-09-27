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

#ifndef SINGULARITY_OPAC_NEUTRINOS_TOPHAT_EMISSIVITY_NEUTRINOS_
#define SINGULARITY_OPAC_NEUTRINOS_TOPHAT_EMISSIVITY_NEUTRINOS_

#include <cassert>
#include <cmath>
#include <cstdio>

#include <ports-of-call/portability.hpp>
#include <singularity-opac/base/opac_error.hpp>
#include <singularity-opac/neutrinos/thermal_distributions_neutrinos.hpp>

namespace singularity {
namespace neutrinos {

// Neutrino tophat emissivity from
// Miller, Ryan, Dolence (2019). arXiv:1903.09273
template <typename ThermalDistribution, typename pc = PhysicalConstantsCGS>
class TophatEmissivity {
 public:
  using PC = pc;

  TophatEmissivity(const Real C, const Real numin, const Real numax)
      : C_(C), numin_(numin), numax_(numax) {}
  TophatEmissivity(const ThermalDistribution &dist, const Real C,
                   const Real numin, const Real numax)
      : dist_(dist), C_(C), numin_(numin), numax_(numax) {}
  TophatEmissivity() = default;
  TophatEmissivity GetOnDevice() { return *this; }
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return 0; }
  PORTABLE_INLINE_FUNCTION
  void PrintParams() const noexcept {
    printf("Tophat emissivity. C, numin, numax = %g, %g, %g\n", C_, numin_,
           numax_);
  }

  PORTABLE_INLINE_FUNCTION
  pc GetPhysicalConstants() const { return pc(); }

  inline void Finalize() noexcept {}

  PORTABLE_INLINE_FUNCTION
  Real AbsorptionCoefficient(const Real rho, const Real temp, const Real Ye,
                             const RadiationType type, const Real nu,
                             Real *lambda = nullptr) const {
    return dist_.AbsorptionCoefficientFromKirkhoff(*this, rho, temp, Ye, type,
                                                   nu, lambda) /
           (4. * M_PI);
  }

  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void
  AbsorptionCoefficient(const Real rho, const Real temp, const Real Ye,
                        const RadiationType type, FrequencyIndexer &nu_bins,
                        DataIndexer &coeffs, const int nbins,
                        Real *lambda = nullptr) const {
    for (int i = 0; i < nbins; ++i) {
      coeffs[i] =
          AbsorptionCoefficient(rho, temp, Ye, type, nu_bins[i], lambda);
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real AngleAveragedAbsorptionCoefficient(const Real rho, const Real temp,
                                          const Real Ye,
                                          const RadiationType type,
                                          const Real nu,
                                          Real *lambda = nullptr) const {
    return dist_.AngleAveragedAbsorptionCoefficientFromKirkhoff(
        *this, rho, temp, Ye, type, nu, lambda);
  }

  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void AngleAveragedAbsorptionCoefficient(
      const Real rho, const Real temp, const Real Ye, const RadiationType type,
      FrequencyIndexer &nu_bins, DataIndexer &coeffs, const int nbins,
      Real *lambda = nullptr) const {
    for (int i = 0; i < nbins; ++i) {
      coeffs[i] = AngleAveragedAbsorptionCoefficient(rho, temp, Ye, type,
                                                     nu_bins[i], lambda);
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real EmissivityPerNuOmega(const Real rho, const Real temp, const Real Ye,
                            const RadiationType type, const Real nu,
                            Real *lambda = nullptr) const {
    if (nu > numin_ && nu < numax_) {
      return C_ * GetYeF(type, Ye);
    } else {
      return 0.;
    }
  }

  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void
  EmissivityPerNuOmega(const Real rho, const Real temp, const Real Ye,
                       const RadiationType type, FrequencyIndexer &nu_bins,
                       DataIndexer &coeffs, const int nbins,
                       Real *lambda = nullptr) const {
    for (int i = 0; i < nbins; ++i) {
      coeffs[i] = EmissivityPerNuOmega(rho, temp, Ye, type, nu_bins[i], lambda);
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real EmissivityPerNu(const Real rho, const Real temp, const Real Ye,
                       const RadiationType type, const Real nu,
                       Real *lambda = nullptr) const {
    return 4 * M_PI * EmissivityPerNuOmega(rho, temp, Ye, type, nu, lambda);
  }

  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void
  EmissivityPerNu(const Real rho, const Real temp, const Real Ye,
                  const RadiationType type, FrequencyIndexer &nu_bins,
                  DataIndexer &coeffs, const int nbins,
                  Real *lambda = nullptr) const {
    for (int i = 0; i < nbins; ++i) {
      coeffs[i] = EmissivityPerNu(rho, temp, Ye, type, nu_bins[i], lambda);
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real Emissivity(const Real rho, const Real temp, const Real Ye,
                  const RadiationType type, Real *lambda = nullptr) const {
    Real Bc = C_ * (numax_ - numin_);
    Real J = 4 * M_PI * Bc * GetYeF(type, Ye);
    return J;
  }

  PORTABLE_INLINE_FUNCTION
  Real NumberEmissivity(const Real rho, const Real temp, const Real Ye,
                        const RadiationType type,
                        Real *lambda = nullptr) const {
    Real Ac = 1 / (pc::h * rho) * C_ * log(numax_ / numin_);
    return 4 * M_PI * rho * Ac * GetYeF(type, Ye);
  }

  PORTABLE_INLINE_FUNCTION
  Real ThermalDistributionOfTNu(const Real temp, const RadiationType type,
                                const Real nu, Real *lambda = nullptr) const {
    return dist_.ThermalDistributionOfTNu(temp, type, nu, lambda);
  }

  PORTABLE_INLINE_FUNCTION
  Real DThermalDistributionOfTNuDT(const Real temp, const RadiationType type,
                                   const Real nu,
                                   Real *lambda = nullptr) const {
    return dist_.DThermalDistributionOfTNuDT(temp, type, nu, lambda);
  }

  PORTABLE_INLINE_FUNCTION
  Real ThermalDistributionOfT(const Real temp, const RadiationType type,
                              Real *lambda = nullptr) const {
    return dist_.ThermalDistributionOfT(temp, type, lambda);
  }

  PORTABLE_INLINE_FUNCTION Real ThermalNumberDistributionOfT(
      const Real temp, const RadiationType type, Real *lambda = nullptr) const {
    return dist_.ThermalNumberDistributionOfT(temp, type, lambda);
  }

  PORTABLE_INLINE_FUNCTION
  Real EnergyDensityFromTemperature(const Real temp, const RadiationType type,
                                    Real *lambda = nullptr) const {
    return dist_.EnergyDensityFromTemperature(temp, type, lambda);
  }

  PORTABLE_INLINE_FUNCTION
  Real TemperatureFromEnergyDensity(const Real er, const RadiationType type,
                                    Real *lambda = nullptr) const {
    return dist_.TemperatureFromEnergyDensity(er, type, lambda);
  }

  PORTABLE_INLINE_FUNCTION
  Real NumberDensityFromTemperature(const Real temp, const RadiationType type,
                                    Real *lambda = nullptr) const {
    return dist_.NumberDensityFromTemperature(temp, type, lambda);
  }

 private:
  PORTABLE_INLINE_FUNCTION
  Real GetYeF(RadiationType type, Real Ye) const {
    if (type == RadiationType::NU_ELECTRON) {
      return 2. * Ye;
    } else if (type == RadiationType::NU_ELECTRON_ANTI) {
      return 1. - 2. * Ye;
    } else {
      return 0.;
    }
  }
  Real C_;
  Real numin_;
  Real numax_;
  ThermalDistribution dist_;
};

} // namespace neutrinos
} // namespace singularity

#endif // SINGULARITY_OPAC_NEUTRINOS_TOPHAT_EMISSIVITY_NEUTRINOS_
