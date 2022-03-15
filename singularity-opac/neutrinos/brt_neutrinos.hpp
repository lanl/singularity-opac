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

#ifndef SINGULARITY_OPAC_NEUTRINOS_BRT_NEUTRINOS_
#define SINGULARITY_OPAC_NEUTRINOS_BRT_NEUTRINOS_

#include <cassert>
#include <cmath>
#include <cstdio>

#include <ports-of-call/portability.hpp>
#include <singularity-opac/base/opac_error.hpp>
#include <singularity-opac/neutrinos/thermal_distributions_neutrinos.hpp>

namespace singularity {
namespace neutrinos {

using pc = PhysicalConstants<CGS>;

// Neutrino electron absorption from Burrows, Reddy, & Thompson 2004
template <typename ThermalDistribution>
class BRTOpacity {
 public:
  BRTOpacity() = default;
  BRTOpacity GetOnDevice() { return *this; }
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return 0; }
  PORTABLE_INLINE_FUNCTION
  void PrintParams() const noexcept {
    printf("Burrows-Reddy-Thompson analytic neutrino opacity.\n");
  }
  inline void Finalize() noexcept {}

  PORTABLE_INLINE_FUNCTION
  Real AbsorptionCoefficient(const Real rho, const Real temp, const Real Ye,
                             const RadiationType type, const Real nu,
                             Real *lambda = nullptr) const {
    return rho / mu_ * GetSigmac(type, nu);
  }
  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void
  AbsorptionCoefficient(const Real rho, const Real temp, const Real Ye,
                        const RadiationType type, FrequencyIndexer &nu_bins,
                        DataIndexer &coeffs, const int nbins,
                        Real *lambda = nullptr) const {
    for (int i = 0; i < nbins; ++i) {
      coeffs[i] = rho / mu_ * GetSigmac(type, nu_bins[i]);
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real AngleAveragedAbsorptionCoefficient(const Real rho, const Real temp,
                                          const Real Ye,
                                          const RadiationType type,
                                          const Real nu,
                                          Real *lambda = nullptr) const {
    return rho / mu_ * GetSigmac(type, nu);
  }

  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void AngleAveragedAbsorptionCoefficient(
      const Real rho, const Real temp, const Real Ye, const RadiationType type,
      FrequencyIndexer &nu_bins, DataIndexer &coeffs, const int nbins,
      Real *lambda = nullptr) const {
    for (int i = 0; i < nbins; ++i) {
      coeffs[i] = rho / mu_ * GetSigmac(type, nu_bins[i]);
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real EmissivityPerNuOmega(const Real rho, const Real temp, const Real Ye,
                            const RadiationType type, const Real nu,
                            Real *lambda = nullptr) const {
    Real Bnu = dist_.ThermalDistributionOfTNu(temp, type, nu, lambda);
    return rho / mu_ * GetSigmac(type, nu) * Bnu;
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
    Real B = dist_.ThermalDistributionOfT(temp, type, lambda);
    return 4 * M_PI * GetAlphac(rho, temp, type);
  }

  PORTABLE_INLINE_FUNCTION
  Real NumberEmissivity(const Real rho, const Real temp, Real Ye,
                        RadiationType type, Real *lambda = nullptr) const {
    return 4 * M_PI * GetNAlphac(rho, temp, type);
  }

  PORTABLE_INLINE_FUNCTION
  Real ThermalDistributionOfTNu(const Real temp, const RadiationType type,
                                const Real nu, Real *lambda = nullptr) const {
    return dist_.ThermalDistributionOfTNu(temp, type, nu, lambda);
  }

  PORTABLE_INLINE_FUNCTION
  Real ThermalDistributionOfT(const Real temp, const RadiationType type,
                              Real *lambda = nullptr) const {
    return dist_.ThermalDistributionOfT(temp, type, lambda);
  }

  PORTABLE_INLINE_FUNCTION Real ThermalNumberDistribution(
      const Real temp, const RadiationType type, Real *lambda = nullptr) const {
    return dist_.ThermalNumberDistribution(temp, type, lambda);
  }

 private:
  Real GetSigmac(const RadiationType type, const Real nu) const {
    if (type != RadiationType::NU_ELECTRON) {
      return 0.;
    }

    return sigma0_ * ((1. + 3. * gA_ * gA_) / 4.) *
           pow((pc::h * nu + Deltanp_) / (pc::me * pc::c * pc::c), 2);
  }

  Real GetAlphac(const Real rho, const Real temp, const RadiationType type) {
    if (type != RadiationType::NU_ELECTRON) {
      return 0.;
    }

    Real retval = (1. + 3. * gA_ * gA_) * pow(pc::kb * temp, 4) * rho * sigma0_;
    retval *= (310. * pow(pc::kb * temp, 2) * pow(M_PI, 6) +
               147. * pow(M_PI, 4) * Deltanp_ * Deltanp_ +
               113400. * pc::kb * temp * Deltanp_ * zeta5_);
    retval /= (5040. * pow(pc::c, 6) * pow(pc::h, 3) * pow(pc::me, 2) * mu_);
    return retval;
  }

  Real GetNAlphac(const Real rho, const Real temp, const RadiationType type) {
    Real retval = (1. + 3. * gA_ * gA_) * pow(pc::kb * temp, 3) * rho * sigma0_;
    retval *= (7. * pc::kb * pow(M_PI, 4) * temp * Deltanp_ +
               90. * pow(Deltanp_, 2) * zeta3_ +
               1350. * pow(pc::kb * temp, 2) * zeta5_);
    retval /= (120. * pow(pc::c, 6) * pow(pc::h, 3) * pow(pc::me, 2) * mu_);
    return retval;
  }

  const Real Fc_ =
      4.543791885043567014e+00; // Fermi coupling constant. Units of erg^-2
  const Real sigma0_ = 4. * Fc_ * Fc_ * pc::c * pc::c * pc::hbar * pc::hbar *
                       pc::me * pc::me * pc::c * pc::c * pc::c * pc::c /
                       M_PI; // Fiducial weak cross section. Units of cm^2
  const Real gA_ = -1.23;
  const Real Deltanp_ = 2.072126995e-6; // erg
  const Real mu_ = pc::mp;
  const Real zeta3_ = 1.2020569031595942853;
  const Real zeta5_ = 1.0369277551433699263;
  ThermalDistribution dist_;
};

} // namespace neutrinos
} // namespace singularity

#endif // SINGULARITY_OPAC_NEUTRINOS_BRT_NEUTRINOS_
