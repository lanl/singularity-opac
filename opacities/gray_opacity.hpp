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

#ifndef OPACITIES_GRAY_EMISSIVITY
#define OPACITIES_GRAY_EMISSIVITY

#include <cassert>
#include <cmath>
#include <cstdio>

#include <opac-utils/opac_error.hpp>
#include <opac-utils/physical_constants.hpp>
#include <ports-of-call/portability.hpp>

#include "thermal_distributions.hpp"

namespace singularity {

template <typename ThermalDistribution, int NSPECIES> class GrayOpacity {
public:
  GrayOpacity(const Real kappa) : kappa_(kappa) {}
  GrayOpacity(const ThermalDistribution &dist, const Real kappa)
      : dist_(dist), kappa_(kappa) {}

  GrayOpacity GetOnDevice() { return *this; }
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return 1; }
  PORTABLE_INLINE_FUNCTION
  void PrintParams() const noexcept {
    printf("Gray opacity. kappa = %g\n", kappa_);
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
    Real Bnu = dist_.ThermalDistributionOfTNu(temp, nu);
    return kappa_ * Bnu;
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
    return kappa_ * dist_.ThermalDistributionOfT(temp);
  }

  PORTABLE_INLINE_FUNCTION
  Real NumberEmissivity(RadiationType type, const Real rho, const Real temp,
                        Real *lambda = nullptr) {
    using namespace constants;
    constexpr Real zeta3 = 1.20206;
    return 12. * pow(cgs::KBOL, 3) * M_PI * NSPECIES * pow(temp, 3) * kappa_ *
           zeta3 / (pow(cgs::CL, 2) * pow(cgs::HPL, 3));
  }

private:
  Real kappa_; // absorption coefficient. Units of 1/cm
  ThermalDistribution dist_;
};

} // namespace singularity

#endif //  OPACITIES_GRAY_EMISSIVITY
