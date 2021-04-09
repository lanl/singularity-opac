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

#include "../utils/ports-of-call/portability.hpp"
#include "../utils/variant/include/mpark/variant.hpp"

using Real = double;

namespace singularity {

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
        [](auto &opac) { return opac_variant<OPACs...>(opac.GetOnDevice()); },
        opac_);
  }

  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept {
    return mpark::visit([](const auto &eos) { return eos.nlambda(); }, eos_);
  }

  template <typename T> PORTABLE_INLINE_FUNCTION bool IsType() const noexcept {
    return mpark::holds_alternative<T>(eos_);
  }

  PORTABLE_INLINE_FUNCTION
  void PrintParams() const noexcept {
    return mpark::visit([](const auto &eos) { return eos.PrintParams(); },
                        eos_);
  }

  inline void Finalize() noexcept {
    return mpark::visit([](auto &eos) { return eos.Finalize(); }, eos_);
  }
}

} // namespace singularity

#endif
