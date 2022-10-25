// ======================================================================
// © 2022. Triad National Security, LLC. All rights reserved.  This
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

#ifndef SINGULARITY_OPAC_PHOTONS_EPBREMSSTRAHLUNG_OPACITY_PHOTONS_
#define SINGULARITY_OPAC_PHOTONS_EPBREMSSTRAHLUNG_OPACITY_PHOTONS_

#include <cassert>
#include <cmath>
#include <cstdio>

#include <ports-of-call/portability.hpp>
#include <singularity-opac/base/opac_error.hpp>
#include <singularity-opac/photons/thermal_distributions_photons.hpp>

namespace singularity {
namespace photons {

template <typename pc = PhysicalConstantsCGS>
class EPBremsstrahlungOpacity {
 public:
  EPBremsstrahlungOpacity() = default;
  EPBremsstrahlungOpacity(const PlanckDistribution<pc> &dist) : dist_(dist) {}

  EPBremsstrahlungOpacity GetOnDevice() { return *this; }
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return 0; }
  PORTABLE_INLINE_FUNCTION
  void PrintParams() const noexcept {
    printf("Electron-proton bremsstrahlung opacity.\n");
  }
  inline void Finalize() noexcept {}

  PORTABLE_INLINE_FUNCTION
  Real AbsorptionCoefficient(const Real rho, const Real temp, const Real nu,
                             Real *lambda = nullptr) const {
    return dist_.AbsorptionCoefficientFromKirkhoff(*this, rho, temp, nu,
                                                   lambda);
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
    return dist_.AngleAveragedAbsorptionCoefficientFromKirkhoff(
        *this, rho, temp, nu, lambda);
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
    const Real thetaE = pc::kb * temp / (pc::me * pc::c * pc::c);
    const Real x = pc::h * nu / (pc::kb * temp);
    const Real ne = rho / mmw_;
    const Real ni = ne;
    return prefac_ / std::sqrt(thetaE) * ne * ni * std::exp(-x) * gff_;
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

    const Real thetaE = pc::kb * temp / (pc::me * pc::c * pc::c);
    const Real ne = rho / mmw_;
    const Real ni = ne;
    return 4. * M_PI * pc::kb * temp / pc::h * prefac_ / std::sqrt(thetaE) *
           ne * ni * gff_;
  }

  PORTABLE_INLINE_FUNCTION
  Real NumberEmissivity(const Real rho, const Real temp,
                        Real *lambda = nullptr) const {
    // Infrared catastrophe
    return std::numeric_limits<Real>::infinity();
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

 private:
  Real mass_ion_;
  PlanckDistribution<pc> dist_;
  static constexpr Real mmw_ =
      (pc::mp + pc::me) /
      (2. * pc::mp); // Neutral fully ionized electron-proton gas
  static constexpr Real gff_ = 1.2;
  static constexpr Real prefac_ = 8. * std::pow(pc::qe, 6) /
                                  (3. * std::pow(pc::me * pc::c * pc::c, 2)) *
                                  std::sqrt(2. * M_PI / 3.);
};

} // namespace photons
} // namespace singularity

#endif // SINGULARITY_OPAC_PHOTONS_EPBREMSSTRAHLUNG_OPACITY_PHOTONS_
