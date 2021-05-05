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

#ifndef OPACITIES_OPAC_VARIANT_HPP_
#define OPACITIES_OPAC_VARIANT_HPP_

#include <opac-utils/opac_error.hpp>
#include <opac-utils/radiation_types.hpp>
#include <ports-of-call/portability.hpp>
#include <variant/include/mpark/variant.hpp>

namespace singularity {

namespace opac_impl {

template <typename... Ts> using opac_variant = mpark::variant<Ts...>;

template <typename... Opacs> class Variant {
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

  // opacity
  PORTABLE_INLINE_FUNCTION
  Real OpacityPerNu(const RadiationType type, const Real rho, const Real temp,
                    const Real nu, Real *lambda = nullptr) {
    mpark::visit(
        [=](const auto &opac) {
          return opac.OpacityPerNu(rho, temp, nu, type, lambda);
        },
        opac_);
  }

  // emissivity with units of energy/time/frequency/volume/angle
  PORTABLE_INLINE_FUNCTION
  Real EmissivityPerNuOmega(const RadiationType type, const Real rho,
                            const Real temp, const Real nu,
                            Real *lambda = nullptr) {
    mpark::visit(
        [=](const auto &opac) {
          return opac.EmissivityPerNuOmega(rho, temp, nu, type, lambda);
        },
        opac_);
  }

  // emissivity integrated over angle
  PORTABLE_INLINE_FUNCTION
  Real EmissivityPerNu(const RadiationType type, const Real rho,
                       const Real temp, const Real nu, Real *lambda = nullptr) {
    mpark::visit(
        [=](const auto &opac) {
          return opac.EmissivityPerNu(rho, temp, nu, type, lambda);
        },
        opac_);
  }

  // emissivity integrated over angle and frequency
  PORTABLE_INLINE_FUNCTION
  Real Emissivity(const RadiationType type, const Real rho, const Real temp,
                  Real *lambda = nullptr) {
    mpark::visit(
        [=](const auto &opac) {
          return opac.Emissivity(rho, temp, type, lambda);
        },
        opac_);
  }

  // Emissivity of packet
  PORTABLE_INLINE_FUNCTION
  Real NumberEmissivity(RadiationType type, const Real rho, const Real temp,
                        Real *lambda = nullptr) {
    mpark::visit(
        [=](const auto &opac) {
          return opac.NumberEmissivity(rho, temp, type, lambda);
        },
        opac_);
  }

  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept {
    return mpark::visit([](const auto &opac) { return opac.nlambda(); }, opac_);
  }

  template <typename T> PORTABLE_INLINE_FUNCTION bool IsType() const noexcept {
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

} // namespace opac_impl
} // namespace singularity

#endif
