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

#ifndef OPACITIES_THERMAL_DISTRIBUTIONS
#define OPACITIES_THERMAL_DISTRIBUTIONS

#include <cmath>

#include <opac-utils/physical_constants.hpp>
#include <opac-utils/radiation_types.hpp>
#include <ports-of-call/portability.hpp>

namespace singularity {

template <int NSPECIES> struct PlanckDistribution {
  PORTABLE_INLINE_FUNCTION
  static Real ThermalDistributionOfTNu(const Real temp, const Real nu) {
    using namespace constants;
    Real x = cgs::HPL * nu / (cgs::KBOL * temp);
    Real Bnu = NSPECIES * (2. * cgs::HPL * nu * nu * nu / (cgs::CL * cgs::CL)) *
               1. / (std::exp(x) + 1.);
    return Bnu;
  }
  PORTABLE_INLINE_FUNCTION
  static Real ThermalDistributionOfT(const Real temp) {
    using namespace constants;
    return 8. * std::pow(M_PI, 5) * std::pow(cgs::KBOL, 4) * NSPECIES *
           std::pow(temp, 4) /
           (15. * std::pow(cgs::CL, 2) * std::pow(cgs::HPL, 3));
  }
};

template <typename Emissivity, typename ThermalDistribution>
Real OpacityFromKirkhoff(Emissivity &J, ThermalDistribution &B,
                         const RadiationType type,
                         const Real rho, const Real temp, const Real nu,
                         Real *lambda=nullptr) {
  Real Bnu = B.ThermalDistributionOfTNu(temp, nu);
  Real jnu = J.EmissivityPerNuOmega(type, rho, temp, nu, lambda);
  return jnu / Bnu;
}

} // namespace singularity

#endif // OPACITIES_THERMAL_DISTRIBUTIONS
