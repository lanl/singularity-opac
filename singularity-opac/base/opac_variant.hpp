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

#ifndef SINGULARITY_OPAC_BASE_OPAC_VARIANT_
#define SINGULARITY_OPAC_BASE_OPAC_VARIANT_

#include <utility>

#include <ports-of-call/portability.hpp>
#include <singularity-opac/base/opac_error.hpp>
#include <singularity-opac/base/radiation_types.hpp>
#include <variant/include/mpark/variant.hpp>

namespace singularity {

namespace opac_impl {

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

  // TODO(JMM): Is this variatic magic too much? Would it be better to be more
  // explicit? This is pretty gross.

  // opacity
  // Signature should be at least
  // rho, temp, nu, lambda
  template <typename... Args>
  PORTABLE_INLINE_FUNCTION auto
  AbsorptionCoefficientPerNu(Args &&...args) const {
    return mpark::visit(
        [=](const auto &opac) {
          return opac.AbsorptionCoefficientPerNu(std::forward<Args>(args)...);
        },
        opac_);
  }

  template <typename FrequencyIndexer, typename DataIndexer, typename... Args>
  PORTABLE_INLINE_FUNCTION void
  AbsorptionCoefficientPerNuBin(const FrequencyIndexer &nu_bins,
                                DataIndexer &coeffs, int nbins,
                                const Real rho, const Real temp, Args&&...args) const {
    mpark::visit(
        [=](const auto &opac) {
          opac.AbsorptionCoefficientPerNuBin(nu_bins, coeffs, nbins,
                                             rho, temp, std::forward<Args>(args)...);
        },
        opac_);
  }

  // emissivity with units of energy/time/frequency/volume/angle
  // signature should be at least rho, temp, nu, lambda
  template <typename... Args>
  PORTABLE_INLINE_FUNCTION auto EmissivityPerNuOmega(Args &&...args) const {
    return mpark::visit(
        [=](const auto &opac) {
          return opac.EmissivityPerNuOmega(std::forward<Args>(args)...);
        },
        opac_);
  }

  // emissivity integrated over angle
  template <typename... Args>
  PORTABLE_INLINE_FUNCTION auto EmissivityPerNu(Args &&...args) const {
    return mpark::visit(
        [=](const auto &opac) {
          return opac.EmissivityPerNu(std::forward<Args>(args)...);
        },
        opac_);
  }

  // emissivity integrated over angle and frequency
  template <typename... Args>
  PORTABLE_INLINE_FUNCTION auto Emissivity(Args &&...args) const {
    return mpark::visit(
        [=](const auto &opac) {
          return opac.Emissivity(std::forward<Args>(args)...);
        },
        opac_);
  }

  // Emissivity of packet
  template <typename... Args>
  PORTABLE_INLINE_FUNCTION auto NumberEmissivity(Args &&...args) const {
    return mpark::visit(
        [=](const auto &opac) {
          return opac.NumberEmissivity(std::forward<Args>(args)...);
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

} // namespace opac_impl
} // namespace singularity

#endif // SINGULARITY_OPAC_BASE_OPAC_VARIANT_
