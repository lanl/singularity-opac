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

#ifndef SINGULARITY_OPAC_PHOTONS_PHOTON_VARIANT_
#define SINGULARITY_OPAC_PHOTONS_PHOTON_VARIANT_

#include <utility>

#include <ports-of-call/portability.hpp>
#include <singularity-opac/base/opac_error.hpp>
#include <singularity-opac/base/radiation_types.hpp>
#include <variant/include/mpark/variant.hpp>

namespace singularity {
namespace photons {
namespace impl {

template <typename... Ts>
using opac_variant = mpark::variant<Ts...>;

template <typename... Opacs>
class Variant {
 private:
  opac_variant<Opacs...> opac_;

 public:
  template <
      typename Choice,
      typename std::enable_if<
          !std::is_same<Variant, typename std::decay<Choice>::type>::value,
          bool>::type = true>
  PORTABLE_FUNCTION Variant(Choice &&choice)
      : opac_(std::forward<Choice>(choice)) {}

  PORTABLE_FUNCTION
  Variant() noexcept = default;

  template <
      typename Choice,
      typename std::enable_if<
          !std::is_same<Variant, typename std::decay<Choice>::type>::value,
          bool>::type = true>
  PORTABLE_FUNCTION Variant &operator=(Choice &&opac) {
    opac_ = std::forward<Choice>(opac);
    return *this;
  }

  template <
      typename Choice,
      typename std::enable_if<
          !std::is_same<Variant, typename std::decay<Choice>::type>::value,
          bool>::type = true>
  Choice get() {
    return mpark::get<Choice>(opac_);
  }

  Variant GetOnDevice() {
    return mpark::visit(
        [](auto &opac) { return opac_variant<Opacs...>(opac.GetOnDevice()); },
        opac_);
  }

  // Directional absorption coefficient with units of 1/length
  // Signature should be at least
  // rho, temp, nu, lambda
  PORTABLE_INLINE_FUNCTION Real
  AbsorptionCoefficient(const Real rho, const Real temp, const Real nu,
                        Real *lambda = nullptr) const {
    return mpark::visit(
        [=](const auto &opac) {
          return opac.AbsorptionCoefficient(rho, temp, nu, lambda);
        },
        opac_);
  }

  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void
  AbsorptionCoefficient(const Real rho, const Real temp,
                        FrequencyIndexer &nu_bins, DataIndexer &coeffs,
                        const int nbins, Real *lambda = nullptr) const {
    mpark::visit(
        [&](const auto &opac) {
          opac.AbsorptionCoefficient(rho, temp, nu_bins, coeffs, nbins, lambda);
        },
        opac_);
  }

  // Angle-averaged absorption coefficient with units of 1/length
  // Signature should be at least
  // rho, temp, nu, lambda
  PORTABLE_INLINE_FUNCTION Real AngleAveragedAbsorptionCoefficient(
      const Real rho, const Real temp, const Real nu,
      Real *lambda = nullptr) const {
    return mpark::visit(
        [=](const auto &opac) {
          return opac.AngleAveragedAbsorptionCoefficient(rho, temp, nu, lambda);
        },
        opac_);
  }

  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void AngleAveragedAbsorptionCoefficient(
      const Real rho, const Real temp, FrequencyIndexer &nu_bins,
      DataIndexer &coeffs, const int nbins, Real *lambda = nullptr) const {
    mpark::visit(
        [&](const auto &opac) {
          opac.AngleAveragedAbsorptionCoefficient(rho, temp, nu_bins, coeffs,
                                                  nbins, lambda);
        },
        opac_);
  }

  // emissivity with units of energy/time/frequency/volume/angle
  // signature should be at least rho, temp, nu, lambda
  PORTABLE_INLINE_FUNCTION Real
  EmissivityPerNuOmega(const Real rho, const Real temp, const Real nu,
                       Real *lambda = nullptr) const {
    return mpark::visit(
        [=](const auto &opac) {
          return opac.EmissivityPerNuOmega(rho, temp, nu, lambda);
        },
        opac_);
  }

  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void
  EmissivityPerNuOmega(const Real rho, const Real temp,
                       FrequencyIndexer &nu_bins, DataIndexer &coeffs,
                       const int nbins, Real *lambda = nullptr) const {
    mpark::visit(
        [&](const auto &opac) {
          return opac.EmissivityPerNuOmega(rho, temp, nu_bins, coeffs, nbins,
                                           lambda);
        },
        opac_);
  }

  // emissivity integrated over angle
  PORTABLE_INLINE_FUNCTION Real EmissivityPerNu(const Real rho, const Real temp,
                                                const Real nu,
                                                Real *lambda = nullptr) const {
    return mpark::visit(
        [=](const auto &opac) {
          return opac.EmissivityPerNu(rho, temp, nu, lambda);
        },
        opac_);
  }

  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void
  EmissivityPerNu(const Real rho, const Real temp, FrequencyIndexer &nu_bins,
                  DataIndexer &coeffs, const int nbins,
                  Real *lambda = nullptr) const {
    mpark::visit(
        [&](const auto &opac) {
          return opac.EmissivityPerNu(rho, temp, nu_bins, coeffs, nbins,
                                      lambda);
        },
        opac_);
  }

  // emissivity integrated over angle and frequency
  PORTABLE_INLINE_FUNCTION Real Emissivity(const Real rho, const Real temp,
                                           Real *lambda = nullptr) const {
    return mpark::visit(
        [=](const auto &opac) { return opac.Emissivity(rho, temp, lambda); },
        opac_);
  }

  // Emissivity of packet
  PORTABLE_INLINE_FUNCTION auto NumberEmissivity(const Real rho,
                                                 const Real temp,
                                                 Real *lambda = nullptr) const {
    return mpark::visit(
        [=](const auto &opac) {
          return opac.NumberEmissivity(rho, temp, lambda);
        },
        opac_);
  }

  // Specific intensity of thermal distribution
  PORTABLE_INLINE_FUNCTION Real ThermalDistributionOfTNu(
      const Real temp, const Real nu, Real *lambda = nullptr) const {
    return mpark::visit(
        [=](const auto &opac) {
          return opac.ThermalDistributionOfTNu(temp, nu, lambda);
        },
        opac_);
  }

  // Integral of thermal distribution over frequency and angle
  PORTABLE_INLINE_FUNCTION Real
  ThermalDistributionOfT(const Real temp, Real *lambda = nullptr) const {
    return mpark::visit(
        [=](const auto &opac) {
          return opac.ThermalDistributionOfT(temp, lambda);
        },
        opac_);
  }

  // Integral of thermal distribution over energy per frequency and angle
  PORTABLE_INLINE_FUNCTION Real
  ThermalNumberDistributionOfT(const Real temp, Real *lambda = nullptr) const {
    return mpark::visit(
        [=](const auto &opac) {
          return opac.ThermalNumberDistributionOfT(temp, lambda);
        },
        opac_);
  }

  // Energy density of thermal distribution
  PORTABLE_INLINE_FUNCTION Real
  EnergyDensityFromTemperature(const Real temp, Real *lambda = nullptr) const {
    return mpark::visit(
        [=](const auto &opac) {
          return opac.EnergyDensityFromTemperature(temp, lambda);
        },
        opac_);
  }

  // Temperature of thermal distribution
  PORTABLE_INLINE_FUNCTION Real
  TemperatureFromEnergyDensity(const Real er, Real *lambda = nullptr) const {
    return mpark::visit(
        [=](const auto &opac) {
          return opac.TemperatureFromEnergyDensity(er, lambda);
        },
        opac_);
  }

  // Number density of thermal distribution
  PORTABLE_INLINE_FUNCTION Real
  NumberDensityFromTemperature(const Real temp, Real *lambda = nullptr) const {
    return mpark::visit(
        [=](const auto &opac) {
          return opac.NumberDensityFromTemperature(temp, lambda);
        },
        opac_);
  }

  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept {
    return mpark::visit([](const auto &opac) { return opac.nlambda(); }, opac_);
  }

  template <typename T>
  PORTABLE_INLINE_FUNCTION bool IsType() const noexcept {
    return mpark::holds_alternative<T>(opac_);
  }

  PORTABLE_INLINE_FUNCTION
  void PrintParams() const noexcept {
    return mpark::visit([](const auto &opac) { return opac.PrintParams(); },
                        opac_);
  }

  inline void Finalize() noexcept {
    return mpark::visit([](auto &opac) { return opac.Finalize(); }, opac_);
  }
};

} // namespace impl
} // namespace photons
} // namespace singularity

#endif // SINGULARITY_OPAC_PHOTONS_PHOTON_VARIANT_
