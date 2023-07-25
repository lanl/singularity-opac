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

#ifndef SINGULARITY_OPAC_PHOTONS_MEAN_PHOTON_S_VARIANT_
#define SINGULARITY_OPAC_PHOTONS_MEAN_PHOTON_S_VARIANT_

#include <utility>

#include <mpark/variant.hpp>
#include <ports-of-call/portability.hpp>
#include <singularity-opac/base/opac_error.hpp>
#include <singularity-opac/base/radiation_types.hpp>
#include <singularity-opac/photons/photon_s_variant.hpp>

namespace singularity {
namespace photons {
namespace impl {

template <typename... SOpacs>
class MeanSVariant {
 private:
  s_opac_variant<SOpacs...> s_opac_;

 public:
  template <
      typename Choice,
      typename std::enable_if<
          !std::is_same<MeanSVariant, typename std::decay<Choice>::type>::value,
          bool>::type = true>
  PORTABLE_FUNCTION MeanSVariant(Choice &&choice)
      : s_opac_(std::forward<Choice>(choice)) {}

  MeanSVariant() = default;

  template <
      typename Choice,
      typename std::enable_if<
          !std::is_same<MeanSVariant, typename std::decay<Choice>::type>::value,
          bool>::type = true>
  PORTABLE_FUNCTION MeanSVariant &operator=(Choice &&s_opac) {
    s_opac_ = std::forward<Choice>(s_opac);
    return *this;
  }

  template <
      typename Choice,
      typename std::enable_if<
          !std::is_same<MeanSVariant, typename std::decay<Choice>::type>::value,
          bool>::type = true>
  Choice get() {
    return mpark::get<Choice>(s_opac_);
  }

  MeanSVariant GetOnDevice() {
    return mpark::visit(
        [](auto &s_opac) {
          return s_opac_variant<SOpacs...>(s_opac.GetOnDevice());
        },
        s_opac_);
  }

  PORTABLE_INLINE_FUNCTION Real
  PlanckMeanTotalScatteringCoefficient(const Real rho, const Real temp) const {
    return mpark::visit(
        [=](const auto &s_opac) {
          return s_opac.PlanckMeanTotalScatteringCoefficient(rho, temp);
        },
        s_opac_);
  }
  PORTABLE_INLINE_FUNCTION Real RosselandMeanTotalScatteringCoefficient(
      const Real rho, const Real temp) const {
    return mpark::visit(
        [=](const auto &s_opac) {
          return s_opac.RosselandMeanTotalScatteringCoefficient(rho, temp);
        },
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

#endif // SINGULARITY_OPAC_PHOTONS_MEAN_PHOTON_S_VARIANT_
