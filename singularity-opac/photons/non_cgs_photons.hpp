// ======================================================================
// © 2021-2024. Triad National Security, LLC. All rights reserved. This
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

#ifndef SINGULARITY_OPAC_PHOTONS_NON_CGS_PHOTONS_
#define SINGULARITY_OPAC_PHOTONS_NON_CGS_PHOTONS_

#include <cassert>
#include <cmath>
#include <cstdio>

#include <ports-of-call/portability.hpp>
#include <singularity-opac/base/opac_error.hpp>
#include <singularity-opac/photons/mean_photon_types.hpp>

namespace singularity {
namespace photons {

template <typename Opac>
class NonCGSUnits {
 public:
  using PC = typename Opac::PC;

  NonCGSUnits() = default;
  NonCGSUnits(Opac &&opac, const Real time_unit, const Real mass_unit,
              const Real length_unit, const Real temp_unit)
      : opac_(std::forward<Opac>(opac)), time_unit_(time_unit),
        mass_unit_(mass_unit), length_unit_(length_unit), temp_unit_(temp_unit),
        rho_unit_(mass_unit_ / (length_unit_ * length_unit_ * length_unit_)),
        freq_unit_(1. / time_unit_),
        inv_emiss_unit_(length_unit_ * time_unit_ * time_unit_ / mass_unit_),
        inv_num_emiss_unit_(length_unit_ * length_unit_ * length_unit_ *
                            time_unit_),
        inv_intensity_unit_(time_unit_ * time_unit_ / mass_unit_),
        inv_energy_dens_unit_(time_unit_ * time_unit_ * length_unit_ /
                              mass_unit_) {}

  auto GetOnDevice() {
    return NonCGSUnits<Opac>(opac_.GetOnDevice(), time_unit_, mass_unit_,
                             length_unit_, temp_unit_);
  }
  inline void Finalize() noexcept { opac_.Finalize(); }

  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return opac_.nlambda(); }

  PORTABLE_INLINE_FUNCTION
  Real AbsorptionCoefficient(const Real rho, const Real temp, const Real nu,
                             Real *lambda = nullptr) const {
    const Real alpha = opac_.AbsorptionCoefficient(
        rho_unit_ * rho, temp_unit_ * temp, nu * freq_unit_, lambda);
    // alpha output in units of 1/cm. Want to convert out of CGS.
    // multiplication by length_unit converts length to cm.
    // division converts length from cm to unit system.
    // thus multiplication converts (1/cm) to unit system.
    return alpha * length_unit_;
  }

  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void
  AbsorptionCoefficient(const Real rho, const Real temp,
                        FrequencyIndexer &nu_bins, DataIndexer &coeffs,
                        const int nbins, Real *lambda = nullptr) const {
    for (int i = 0; i < nbins; ++i) {
      nu_bins[i] *= freq_unit_;
    }
    opac_.AbsorptionCoefficient(rho * rho_unit_, temp * temp_unit_, nu_bins,
                                coeffs, nbins, lambda);
    for (int i = 0; i < nbins; ++i) {
      nu_bins[i] *= time_unit_;
      coeffs[i] *= length_unit_;
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real AngleAveragedAbsorptionCoefficient(const Real rho, const Real temp,
                                          const Real nu,
                                          Real *lambda = nullptr) const {
    const Real alpha = opac_.AngleAveragedAbsorptionCoefficient(
        rho * rho_unit_, temp * temp_unit_, nu * freq_unit_, lambda);
    return alpha * length_unit_;
  }

  // TODO(JMM): Doing the frequency conversion here is SUPER gross.
  // Is there no better way? Should we just not allow modified units?
  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void AngleAveragedAbsorptionCoefficient(
      const Real rho, const Real temp, FrequencyIndexer &nu_bins,
      DataIndexer &coeffs, const int nbins, Real *lambda = nullptr) const {
    for (int i = 0; i < nbins; ++i) {
      nu_bins[i] *= freq_unit_;
    }
    opac_.AngleAveragedAbsorptionCoefficient(rho * rho_unit_, temp * temp_unit_,
                                             nu_bins, coeffs, nbins, lambda);
    for (int i = 0; i < nbins; ++i) {
      nu_bins[i] *= time_unit_;
      coeffs[i] *= length_unit_;
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real EmissivityPerNuOmega(const Real rho, const Real temp, const Real nu,
                            Real *lambda = nullptr) const {
    const Real jnu = opac_.EmissivityPerNuOmega(
        rho * rho_unit_, temp * temp_unit_, nu * freq_unit_, lambda);
    return jnu * inv_emiss_unit_;
  }

  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void
  EmissivityPerNuOmega(const Real rho, const Real temp,
                       FrequencyIndexer &nu_bins, DataIndexer &coeffs,
                       const int nbins, Real *lambda = nullptr) const {
    for (int i = 0; i < nbins; ++i) {
      nu_bins[i] *= freq_unit_;
    }
    opac_.EmissivityPerNuOmega(rho * rho_unit_, temp * temp_unit_, nu_bins,
                               coeffs, nbins, lambda);
    for (int i = 0; i < nbins; ++i) {
      nu_bins[i] *= time_unit_;
      coeffs[i] *= inv_emiss_unit_;
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real EmissivityPerNu(const Real rho, const Real temp, const Real nu,
                       Real *lambda = nullptr) const {
    Real Jnu = opac_.EmissivityPerNu(rho * rho_unit_, temp * temp_unit_,
                                     nu * freq_unit_, lambda);
    return Jnu * inv_emiss_unit_; // solid angle is unitless
  }

  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void
  EmissivityPerNu(const Real rho, const Real temp, FrequencyIndexer &nu_bins,
                  DataIndexer &coeffs, const int nbins,
                  Real *lambda = nullptr) const {
    for (int i = 0; i < nbins; ++i) {
      nu_bins[i] *= freq_unit_;
    }
    opac_.EmissivityPerNu(rho * rho_unit_, temp * temp_unit_, nu_bins, coeffs,
                          nbins, lambda);
    for (int i = 0; i < nbins; ++i) {
      nu_bins[i] *= time_unit_;
      coeffs[i] *= inv_emiss_unit_;
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real Emissivity(const Real rho, const Real temp,
                  Real *lambda = nullptr) const {
    const Real J = opac_.Emissivity(rho * rho_unit_, temp * temp_unit_, lambda);
    // Jnu integrated over frequency, but divide by frequency to get out of cgs
    return J * inv_emiss_unit_ * time_unit_;
  }

  PORTABLE_INLINE_FUNCTION
  Real NumberEmissivity(const Real rho, const Real temp,
                        Real *lambda = nullptr) const {
    Real JoH =
        opac_.NumberEmissivity(rho * rho_unit_, temp * temp_unit_, lambda);
    return JoH * inv_num_emiss_unit_;
  }

  PORTABLE_INLINE_FUNCTION
  Real ThermalDistributionOfTNu(const Real temp, const Real nu,
                                Real *lambda = nullptr) const {
    Real BoH = opac_.ThermalDistributionOfTNu(temp * temp_unit_,
                                              nu * freq_unit_, lambda);
    return BoH * inv_intensity_unit_;
  }

  PORTABLE_INLINE_FUNCTION
  Real DThermalDistributionOfTNuDT(const Real temp, const Real nu,
                                   Real *lambda = nullptr) const {
    Real dBdToH = opac_.DThermalDistributionOfTNuDT(temp, nu, lambda);
    return dBdToH * inv_intensity_unit_ * temp_unit_;
  }

  PORTABLE_INLINE_FUNCTION
  Real ThermalDistributionOfT(const Real temp, Real *lambda = nullptr) const {
    Real BoH = opac_.ThermalDistributionOfT(temp * temp_unit_, lambda);
    return BoH * inv_energy_dens_unit_ * time_unit_ / length_unit_;
  }

  PORTABLE_INLINE_FUNCTION
  Real ThermalNumberDistributionOfT(const Real temp,
                                    Real *lambda = nullptr) const {
    Real NoH = opac_.ThermalNumberDistributionOfT(temp * temp_unit_, lambda);
    return NoH * mass_unit_ / rho_unit_ * time_unit_ / length_unit_;
  }

  PORTABLE_INLINE_FUNCTION
  Real EnergyDensityFromTemperature(const Real temp,
                                    Real *lambda = nullptr) const {
    Real BoH = opac_.EnergyDensityFromTemperature(temp * temp_unit_, lambda);
    return BoH * inv_energy_dens_unit_;
  }

  PORTABLE_INLINE_FUNCTION
  Real TemperatureFromEnergyDensity(const Real er,
                                    Real *lambda = nullptr) const {
    Real BoH =
        opac_.TemperatureFromEnergyDensity(er / inv_energy_dens_unit_, lambda);
    return BoH / temp_unit_;
  }

  PORTABLE_INLINE_FUNCTION
  Real NumberDensityFromTemperature(const Real temp,
                                    Real *lambda = nullptr) const {
    Real NoH = opac_.NumberDensityFromTemperature(temp * temp_unit_, lambda);
    return NoH * mass_unit_ / rho_unit_;
  }

  PORTABLE_INLINE_FUNCTION RuntimePhysicalConstants
  GetRuntimePhysicalConstants() const {
    return RuntimePhysicalConstants(PC(), time_unit_, mass_unit_, length_unit_,
                                    temp_unit_);
  }

 private:
  Opac opac_;
  Real time_unit_, mass_unit_, length_unit_, temp_unit_;
  Real rho_unit_, freq_unit_, inv_emiss_unit_, inv_num_emiss_unit_;
  Real inv_intensity_unit_, inv_energy_dens_unit_;
};

template <typename MeanOpac>
class MeanNonCGSUnits {
 public:
  using PC = typename MeanOpac::PC;

  MeanNonCGSUnits() = default;
  MeanNonCGSUnits(MeanOpac &&mean_opac, const Real time_unit,
                  const Real mass_unit, const Real length_unit,
                  const Real temp_unit)
      : mean_opac_(std::forward<MeanOpac>(mean_opac)), time_unit_(time_unit),
        mass_unit_(mass_unit), length_unit_(length_unit), temp_unit_(temp_unit),
        rho_unit_(mass_unit_ / (length_unit_ * length_unit_ * length_unit_)),
        inv_emiss_unit_(length_unit_ * time_unit_ * time_unit_ * time_unit_ /
          mass_unit_) {}

  auto GetOnDevice() {
    return MeanNonCGSUnits<MeanOpac>(mean_opac_.GetOnDevice(), time_unit_,
                                     mass_unit_, length_unit_, temp_unit_);
  }
  inline void Finalize() noexcept { mean_opac_.Finalize(); }

  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return mean_opac_.nlambda(); }

#ifdef SPINER_USE_HDF
  void Save(const std::string &filename) const {
    return mean_opac_.Save(filename);
  }
#endif

  PORTABLE_INLINE_FUNCTION
  Real PlanckMeanAbsorptionCoefficient(const Real rho, const Real temp) const {
    const Real alpha = mean_opac_.PlanckMeanAbsorptionCoefficient(
        rho_unit_ * rho, temp_unit_ * temp);
    // alpha output in units of 1/cm. Want to convert out of CGS.
    // multiplication by length_unit converts length to cm.
    // division converts length from cm to unit system.
    // thus multiplication converts (1/cm) to unit system.
    return alpha * length_unit_;
  }

  PORTABLE_INLINE_FUNCTION
  Real RosselandMeanAbsorptionCoefficient(const Real rho,
                                          const Real temp) const {
    const Real alpha = mean_opac_.RosselandMeanAbsorptionCoefficient(
        rho_unit_ * rho, temp_unit_ * temp);
    // alpha output in units of 1/cm. Want to convert out of CGS.
    // multiplication by length_unit converts length to cm.
    // division converts length from cm to unit system.
    // thus multiplication converts (1/cm) to unit system.
    return alpha * length_unit_;
  }

  PORTABLE_INLINE_FUNCTION
  Real AbsorptionCoefficient(const Real rho, const Real temp,
                            const int gmode = Rosseland) const {
    const Real alpha = mean_opac_.AbsorptionCoefficient(rho, temp, gmode);
    return alpha * length_unit_;
  }

  PORTABLE_INLINE_FUNCTION
  Real Emissivity(const Real rho, const Real temp,
                  const int gmode = Rosseland,
                  Real *lambda = nullptr) const {
    const Real J = mean_opac_.Emissivity(rho, temp, gmode);
    return J * inv_emiss_unit_;
  }

  PORTABLE_INLINE_FUNCTION RuntimePhysicalConstants
  GetRuntimePhysicalConstants() const {
    return RuntimePhysicalConstants(PhysicalConstantsCGS(), time_unit_,
                                    mass_unit_, length_unit_, temp_unit_);
  }

 private:
  MeanOpac mean_opac_;
  Real time_unit_, mass_unit_, length_unit_, temp_unit_;
  Real rho_unit_, inv_emiss_unit_;
};

} // namespace photons
} // namespace singularity

#endif // SINGULARITY_OPAC_PHOTONS_NON_CGS_PHOTONS_
