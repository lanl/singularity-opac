// ======================================================================
// © 2022-2024. Triad National Security, LLC. All rights reserved.  This
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

#ifndef SINGULARITY_OPAC_PHOTONS_MEAN_PHOTON_VARIANT_
#define SINGULARITY_OPAC_PHOTONS_MEAN_PHOTON_VARIANT_

#include <utility>

#include <ports-of-call/portability.hpp>
#include <singularity-opac/base/opac_error.hpp>
#include <singularity-opac/base/radiation_types.hpp>
#include <singularity-opac/photons/photon_variant.hpp>
#include <variant/include/mpark/variant.hpp>

namespace singularity {
namespace photons {
namespace impl {

template <typename... Opacs>
class MeanVariant {
 private:
  opac_variant<Opacs...> opac_;

 public:
  template <
      typename Choice,
      typename std::enable_if<
          !std::is_same<MeanVariant, typename std::decay<Choice>::type>::value,
          bool>::type = true>
  PORTABLE_FUNCTION MeanVariant(Choice &&choice)
      : opac_(std::forward<Choice>(choice)) {}

  MeanVariant() = default;

  template <
      typename Choice,
      typename std::enable_if<
          !std::is_same<MeanVariant, typename std::decay<Choice>::type>::value,
          bool>::type = true>
  PORTABLE_FUNCTION MeanVariant &operator=(Choice &&opac) {
    opac_ = std::forward<Choice>(opac);
    return *this;
  }

  template <
      typename Choice,
      typename std::enable_if<
          !std::is_same<MeanVariant, typename std::decay<Choice>::type>::value,
          bool>::type = true>
  Choice get() {
    return mpark::get<Choice>(opac_);
  }

  MeanVariant GetOnDevice() {
    return mpark::visit(
        [](auto &opac) { return opac_variant<Opacs...>(opac.GetOnDevice()); },
        opac_);
  }

  PORTABLE_INLINE_FUNCTION Real
  PlanckMeanAbsorptionCoefficient(const Real rho, const Real temp) const {
    return mpark::visit(
        [=](const auto &opac) {
          return opac.PlanckMeanAbsorptionCoefficient(rho, temp);
        },
        opac_);
  }
  PORTABLE_INLINE_FUNCTION Real
  RosselandMeanAbsorptionCoefficient(const Real rho, const Real temp) const {
    return mpark::visit(
        [=](const auto &opac) {
          return opac.RosselandMeanAbsorptionCoefficient(rho, temp);
        },
        opac_);
  }

  inline void Finalize() noexcept {
    return mpark::visit([](auto &opac) { return opac.Finalize(); }, opac_);
  }

  PORTABLE_INLINE_FUNCTION RuntimePhysicalConstants
  GetRuntimePhysicalConstants() const {
    return mpark::visit(
        [=](const auto &opac) {
          using PC = typename std::decay_t<decltype(opac)>::PC;
          return singularity::GetRuntimePhysicalConstants(PC());
        },
        opac_);
  }
};

} // namespace impl
} // namespace photons
} // namespace singularity

#endif // SINGULARITY_OPAC_PHOTONS_MEAN_PHOTON_VARIANT_
