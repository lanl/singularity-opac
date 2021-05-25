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

#ifndef SINGULARITY_OPAC_PHOTONS_GRAY_OPACITY_PHOTONS_
#define SINGULARITY_OPAC_PHOTONS_GRAY_OPACITY_PHOTONS_

#include <cassert>
#include <cmath>
#include <cstdio>

#include <ports-of-call/portability.hpp>
#include <singularity-opac/base/opac_error.hpp>
#include <singularity-opac/photons/thermal_distributions_photons.hpp>

namespace singularity {
namespace photons {

class GrayOpacity {
 public:
  GrayOpacity(const Real kappa) : kappa_(kappa) {}
  GrayOpacity(const PlanckDistribution &dist, const Real kappa)
      : dist_(dist), kappa_(kappa) {}

  GrayOpacity GetOnDevice() { return *this; }
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return 0; }
  PORTABLE_INLINE_FUNCTION
  void PrintParams() const noexcept {
    printf("Gray opacity. kappa = %g\n", kappa_);
  }
  inline void Finalize() noexcept {}

  PORTABLE_INLINE_FUNCTION
  Real AbsorptionCoefficientPerNu(const Real rho, const Real temp,
                                  const Real nu, Real *lambda = nullptr) const {
    return dist_.AbsorptionCoefficientFromKirkhoff(*this, rho, temp, nu,
                                                   lambda);
  }

  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void AbsorptionCoefficientPerNu(
      const Real rho, const Real temp, const FrequencyIndexer &nu_bins,
      DataIndexer &coeffs, const int nbins, Real *lambda = nullptr) const {
    for (int i = 0; i < nbins; ++i) {
      coeffs[i] = AbsorptionCoefficientPerNu(rho, temp, nu_bins[i], lambda);
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real EmissivityPerNuOmega(const Real rho, const Real temp, const Real nu,
                            Real *lambda = nullptr) const {
    Real Bnu = dist_.ThermalDistributionOfTNu(temp, nu, lambda);
    return rho * kappa_ * Bnu;
  }

  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void
  EmissivityPerNuOmega(const Real rho, const Real temp,
                       const FrequencyIndexer &nu_bins, DataIndexer &coeffs,
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
  EmissivityPerNu(const Real rho, const Real temp,
                  const FrequencyIndexer &nu_bins, DataIndexer &coeffs,
                  const int nbins, Real *lambda = nullptr) const {
    for (int i = 0; i < nbins; ++i) {
      coeffs[i] = EmissivityPerNu(rho, temp, nu_bins[i], lambda);
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real Emissivity(const Real rho, const Real temp,
                  Real *lambda = nullptr) const {
    Real B = dist_.ThermalDistributionOfT(temp, lambda);
    return rho * kappa_ * B;
  }

  PORTABLE_INLINE_FUNCTION
  Real NumberEmissivity(const Real rho, const Real temp,
                        Real *lambda = nullptr) const {
    return kappa_ * dist_.ThermalNumberDistribution(temp, lambda);
  }

 private:
  Real kappa_; // Opacity. Units of cm^2/g
  PlanckDistribution dist_;
};

} // namespace photons
} // namespace singularity

#endif // SINGULARITY_OPAC_PHOTONS_GRAY_OPACITY_PHOTONS_
