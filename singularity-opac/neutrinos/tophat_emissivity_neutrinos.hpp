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

using pc = PhysicalConstants<CGS>;

// Neutrino tophat emissivity from
// Miller, Ryan, Dolence (2019). arXiv:1903.09273
template <typename ThermalDistribution>
class TophatEmissivity {
 public:
  TophatEmissivity(const Real C, const Real numin, const Real numax)
      : C_(C), numin_(numin), numax_(numax) {}
  TophatEmissivity(const ThermalDistribution &dist, const Real C,
                   const Real numin, const Real numax)
      : dist_(dist), C_(C), numin_(numin), numax_(numax) {}
  TophatEmissivity GetOnDevice() { return *this; }
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return 0; }
  PORTABLE_INLINE_FUNCTION
  void PrintParams() const noexcept {
    printf("Tophat emissivity. C, numin, numax = %g, %g, %g\n", C_, numin_,
           numax_);
  }
  inline void Finalize() noexcept {}

  // TODO(JMM): Does this make sense for the tophat?
  PORTABLE_INLINE_FUNCTION
  Real AbsorptionCoefficientPerNu(const RadiationType type, const Real rho,
                                  const Real temp, const Real Ye, const Real nu,
                                  Real *lambda = nullptr) const {
    return dist_.AbsorptionCoefficientFromKirkhoff(*this, type, rho, temp, Ye,
                                                   nu, lambda);
  }

  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void AbsorptionCoefficientPerNu(
      const RadiationType type, const Real rho, const Real temp, const Real Ye,
      const FrequencyIndexer &nu_bins, DataIndexer &coeffs, const int nbins,
      Real *lambda = nullptr) const {
    for (int i = 0; i < nbins; ++i) {
      coeffs[i] =
          AbsorptionCoefficientPerNu(type, rho, temp, Ye, nu_bins[i], lambda);
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real EmissivityPerNuOmega(const RadiationType type, const Real rho,
                            const Real temp, const Real Ye, const Real nu,
                            Real *lambda = nullptr) const {
    if (nu > numin_ && nu < numax_) {
      return C_ * GetYeF(type, Ye) / (4. * M_PI);
    } else {
      return 0.;
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real EmissivityPerNu(const RadiationType type, const Real rho,
                       const Real temp, const Real Ye, const Real nu,
                       Real *lambda = nullptr) const {
    return 4 * M_PI * EmissivityPerNuOmega(type, rho, temp, Ye, nu, lambda);
  }

  PORTABLE_INLINE_FUNCTION
  Real Emissivity(const RadiationType type, const Real rho, const Real temp,
                  const Real Ye, Real *lambda = nullptr) const {
    Real Bc = C_ * (numax_ - numin_);
    Real J = Bc * GetYeF(type, Ye);
    return J;
  }

  PORTABLE_INLINE_FUNCTION
  Real NumberEmissivity(const RadiationType type, const Real rho,
                        const Real temp, const Real Ye,
                        Real *lambda = nullptr) const {
    Real Ac = 1 / (pc::h * rho) * C_ * log(numax_ / numin_);
    return rho * Ac * GetYeF(type, Ye);
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
