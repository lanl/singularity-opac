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

#ifndef SINGULARITY_OPAC_PHOTONS_PHOTON_S_VARIANT_
#define SINGULARITY_OPAC_PHOTONS_PHOTON_S_VARIANT_

#include <utility>

#include <ports-of-call/portability.hpp>
#include <singularity-opac/base/opac_error.hpp>
#include <singularity-opac/base/radiation_types.hpp>
#include <variant/include/mpark/variant.hpp>

namespace singularity {
namespace photons {
namespace impl {

template <typename... Ts>
using s_opac_variant = mpark::variant<Ts...>;

template <typename... S_Opacs>
class S_Variant {
 private:
  s_opac_variant<S_Opacs...> s_opac_;

 public:
  template <
      typename Choice,
      typename std::enable_if<
          !std::is_same<S_Variant, typename std::decay<Choice>::type>::value,
          bool>::type = true>
  PORTABLE_FUNCTION S_Variant(Choice &&choice)
      : s_opac_(std::forward<Choice>(choice)) {}

  S_Variant() = default;

  template <
      typename Choice,
      typename std::enable_if<
          !std::is_same<S_Variant, typename std::decay<Choice>::type>::value,
          bool>::type = true>
  PORTABLE_FUNCTION S_Variant &operator=(Choice &&s_opac) {
    s_opac_ = std::forward<Choice>(s_opac);
    return *this;
  }

  template <
      typename Choice,
      typename std::enable_if<
          !std::is_same<S_Variant, typename std::decay<Choice>::type>::value,
          bool>::type = true>
  Choice get() {
    return mpark::get<Choice>(s_opac_);
  }

  S_Variant GetOnDevice() {
    return mpark::visit(
        [](auto &s_opac) {
          return s_opac_variant<S_Opacs...>(s_opac.GetOnDevice());
        },
        s_opac_);
  }

  // Total cross section with units of length^2
  // Signature should be at least
  // rho, temp, Ye, type, nu, lambda
  PORTABLE_INLINE_FUNCTION Real
  TotalCrossSection(const Real rho, const Real temp, const Real nu,
                    Real *lambda = nullptr) const {
    return mpark::visit(
        [=](const auto &s_opac) {
          return s_opac.TotalCrossSection(rho, temp, nu, lambda);
        },
        s_opac_);
  }

  // Differential cross section with units of length^2/steradian
  // Signature should be at least
  // rho, temp, Ye, type, nu, mu, lambda
  PORTABLE_INLINE_FUNCTION Real
  DifferentialCrossSection(const Real rho, const Real temp, const Real nu,
                           const Real mu, Real *lambda = nullptr) const {
    return mpark::visit(
        [=](const auto &s_opac) {
          return s_opac.DifferentialCrossSection(rho, temp, nu, mu, lambda);
        },
        s_opac_);
  }

  // Total scattering coefficient with units of 1/length
  // Signature should be at least
  // rho, temp, Ye, type, nu, lambda
  PORTABLE_INLINE_FUNCTION Real
  TotalScatteringCoefficient(const Real rho, const Real temp, const Real nu,
                             Real *lambda = nullptr) const {
    return mpark::visit(
        [=](const auto &s_opac) {
          return s_opac.TotalScatteringCoefficient(rho, temp, nu, lambda);
        },
        s_opac_);
  }

  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept {
    return mpark::visit([](const auto &s_opac) { return s_opac.nlambda(); },
                        s_opac_);
  }

  template <typename T>
  PORTABLE_INLINE_FUNCTION bool IsType() const noexcept {
    return mpark::holds_alternative<T>(s_opac_);
  }

  PORTABLE_INLINE_FUNCTION
  void PrintParams() const noexcept {
    return mpark::visit([](const auto &s_opac) { return s_opac.PrintParams(); },
                        s_opac_);
  }

  inline void Finalize() noexcept {
    return mpark::visit([](auto &s_opac) { return s_opac.Finalize(); },
                        s_opac_);
  }
};

} // namespace impl
} // namespace photons
} // namespace singularity

#endif // SINGUALRITY_OPAC_PHOTONS_PHOTON_S_VARIANT_
