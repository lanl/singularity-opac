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

#ifndef SINGULARITY_OPAC_NEUTRINOS_THERMAL_DISTRIBUTIONS_NEUTRINOS_
#define SINGULARITY_OPAC_NEUTRINOS_THERMAL_DISTRIBUTIONS_NEUTRINOS_

#include <cmath>

#include <ports-of-call/portability.hpp>
#include <singularity-opac/base/radiation_types.hpp>
#include <singularity-opac/constants/constants.hpp>

namespace singularity {
namespace neutrinos {

using pc = PhysicalConstants<CGS>;

#define EPS (10.0 * std::numeric_limits<Real>::epsilon())

template <int NSPECIES>
struct FermiDiracDistributionNoMu {
  PORTABLE_INLINE_FUNCTION
  Real ThermalDistributionOfTNu(const Real temp, const RadiationType type,
                                const Real nu, Real *lambda = nullptr) const {
    Real x = pc::h * nu / (pc::kb * temp);
    Real Bnu = NSPECIES * (2. * pc::h * nu * nu * nu / (pc::c * pc::c)) * 1. /
               (std::exp(x) + 1.);
    return Bnu;
  }
  PORTABLE_INLINE_FUNCTION
  Real ThermalDistributionOfT(const Real temp, const RadiationType type,
                              Real *lambda = nullptr) const {
    return 7. * std::pow(M_PI, 5) * std::pow(pc::kb, 4) * NSPECIES *
           std::pow(temp, 4) / (15. * std::pow(pc::c, 3) * std::pow(pc::h, 3));
  }
  PORTABLE_INLINE_FUNCTION
  Real ThermalNumberDistribution(const Real temp, const RadiationType type,
                                 Real *lambda = nullptr) const {
    constexpr Real zeta3 = 1.20206;
    return 16. * pow(pc::kb, 3) * M_PI * NSPECIES * pow(temp, 3) * zeta3 /
           (pow(pc::c, 3) * pow(pc::h, 3));
  }
  template <typename Emissivity>
  PORTABLE_INLINE_FUNCTION Real AbsorptionCoefficientFromKirkhoff(
      const Emissivity &J, const Real rho, const Real temp, const Real Ye,
      const RadiationType type, const Real nu, Real *lambda = nullptr) const {
    const Real Bnu =
        std::max(ThermalDistributionOfTNu(temp, type, nu, lambda), EPS);
    const Real jnu =
        std::max(J.EmissivityPerNuOmega(rho, temp, Ye, type, nu, lambda), EPS);
    return jnu / Bnu;
  }
  template <typename Emissivity>
  PORTABLE_INLINE_FUNCTION Real AngleAveragedAbsorptionCoefficientFromKirkhoff(
      const Emissivity &J, const Real rho, const Real temp, const Real Ye,
      const RadiationType type, const Real nu, Real *lambda = nullptr) const {
    const Real Bnu =
        std::max(ThermalDistributionOfTNu(temp, type, nu, lambda), EPS);
    const Real jnu =
        std::max(J.EmissivityPerNu(rho, temp, Ye, type, nu, lambda), EPS) /
        (4. * M_PI);
    return jnu / Bnu;
  }
};

#undef EPS

} // namespace neutrinos
} // namespace singularity

#endif // SINGULARITY_OPAC_NEUTRINOS_THERMAL_DISTRIBUTIONS_NEUTRINOS_
