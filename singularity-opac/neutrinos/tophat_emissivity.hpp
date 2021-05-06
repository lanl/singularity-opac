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

#ifndef OPACITIES_TOPHAT_EMISSIVITY
#define OPACITIES_TOPHAT_EMISSIVITY

#include <cassert>
#include <cmath>
#include <cstdio>

#include <singularity-opac/base/opac_error.hpp>
#include <singularity-opac/base/physical_constants.hpp>
#include <ports-of-call/portability.hpp>

#include "thermal_distributions.hpp"

namespace singularity {

// Neutrino tophat emissivity from
// Miller, Ryan, Dolence (2019). arXiv:1903.09273
template <typename ThermalDistribution> class TophatEmissivity {
public:
  TophatEmissivity(const Real C, const Real numin, const Real numax)
      : C_(C), numin_(numin), numax_(numax) {}
  TophatEmissivity(const ThermalDistribution &dist, const Real C,
                   const Real numin, const Real numax)
      : dist_(dist), C_(C), numin_(numin), numax_(numax) {}
  TophatEmissivity GetOnDevice() { return *this; }
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return 1; }
  PORTABLE_INLINE_FUNCTION
  void PrintParams() const noexcept {
    printf("Tophat emissivity. C, numin, numax = %g, %g, %g\n", C_, numin_,
           numax_);
  }
  inline void Finalize() noexcept {}

  // TODO(JMM): Does this make sense for the tophat?
  PORTABLE_INLINE_FUNCTION
  Real OpacityPerNu(const RadiationType type, const Real rho, const Real temp,
                    const Real nu, Real *lambda = nullptr) {
    return OpacityFromKirkhoff(*this, dist_);
  }

  PORTABLE_INLINE_FUNCTION
  Real EmissivityPerNuOmega(const RadiationType type, const Real rho,
                            const Real temp, const Real nu,
                            Real *lambda = nullptr) {
    assert(type == RadiationType::NU_ELECTRON ||
           type == RadiationType::NU_ELECTRON_ANTI ||
           type == RadiationType::NU_HEAVY);
    Real Ye = lambda[0];
    if (nu > numin_ && nu < numax_) {
      return C_ * GetYeF(Ye, type) / (4. * M_PI);
    } else {
      return 0.;
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real EmissivityPerNu(const RadiationType type, const Real rho,
                       const Real temp, const Real nu,
                       const Real *lambda = nullptr) {
    return 4 * M_PI * EmissivityPerNuOmega(type, rho, temp, nu, lambda);
  }

  PORTABLE_INLINE_FUNCTION
  Real Emissivity(const RadiationType type, const Real rho, const Real temp,
                  Real *lambda = nullptr) {
    assert(type == RadiationType::NU_ELECTRON ||
           type == RadiationType::NU_ELECTRON_ANTI ||
           type == RadiationType::NU_HEAVY);
    Real Ye = lambda[0];
    Real Bc = C_ * (numax_ - numin_);
    Real J = Bc * GetYeF(Ye, type);
    return J;
  }

  PORTABLE_INLINE_FUNCTION
  Real NumberEmissivity(RadiationType type, const Real rho, const Real temp,
                        Real *lambda = nullptr) {
    using namespace constants;
    assert(type == RadiationType::NU_ELECTRON ||
           type == RadiationType::NU_ELECTRON_ANTI ||
           type == RadiationType::NU_HEAVY);
    Real Ye = lambda[0];
    Real Ac = 1 / (cgs::HPL * rho) * C_ * log(numax_ / numin_);
    return rho * Ac * GetYeF(Ye, type);
  }

private:
  PORTABLE_INLINE_FUNCTION
  Real GetYeF(Real Ye, RadiationType type) {
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

} // namespace singularity

#endif // OPACITIES_TOPHAT_EMISSIVITY
