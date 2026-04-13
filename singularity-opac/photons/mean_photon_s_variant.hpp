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

#include <ports-of-call/portability.hpp>
#include <ports-of-call/variant.hpp>
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
    return PortsOfCall::get<Choice>(s_opac_);
  }

  MeanSVariant GetOnDevice() {
    return PortsOfCall::visit(
        [](auto &s_opac) {
          return s_opac_variant<SOpacs...>(s_opac.GetOnDevice());
        },
        s_opac_);
  }

  PORTABLE_INLINE_FUNCTION Real
  PlanckMeanTotalScatteringCoefficient(const Real rho, const Real temp) const {
    return PortsOfCall::visit(
        [=](const auto &s_opac) {
          return s_opac.PlanckMeanTotalScatteringCoefficient(rho, temp);
        },
        s_opac_);
  }
  PORTABLE_INLINE_FUNCTION Real RosselandMeanTotalScatteringCoefficient(
      const Real rho, const Real temp) const {
    return PortsOfCall::visit(
        [=](const auto &s_opac) {
          return s_opac.RosselandMeanTotalScatteringCoefficient(rho, temp);
        },
        s_opac_);
  }

  PORTABLE_INLINE_FUNCTION
  int ngroups() const noexcept {
    return PortsOfCall::visit([](const auto &s_opac) { return s_opac.ngroups(); },
                              s_opac_);
  }

  PORTABLE_INLINE_FUNCTION
  bool HasGroupBounds() const noexcept {
    return PortsOfCall::visit(
        [](const auto &s_opac) { return s_opac.HasGroupBounds(); }, s_opac_);
  }

  PORTABLE_INLINE_FUNCTION
  Real PlanckGroupScatteringCoefficient(const Real rho, const Real temp,
                                        const int group) const {
    return ScatteringCoefficient(rho, temp, group, Planck);
  }

  PORTABLE_INLINE_FUNCTION
  Real RosselandGroupScatteringCoefficient(const Real rho, const Real temp,
                                           const int group) const {
    return ScatteringCoefficient(rho, temp, group, Rosseland);
  }

  PORTABLE_INLINE_FUNCTION
  Real ScatteringCoefficient(const Real rho, const Real temp, const int group,
                             const int gmode = Rosseland) const {
    return PortsOfCall::visit(
        [=](const auto &s_opac) {
          return s_opac.ScatteringCoefficient(rho, temp, group, gmode);
        },
        s_opac_);
  }

  PORTABLE_INLINE_FUNCTION
  int GroupOfNu(const Real nu) const {
    return PortsOfCall::visit(
        [=](const auto &s_opac) { return s_opac.GroupOfNu(nu); }, s_opac_);
  }

  PORTABLE_INLINE_FUNCTION
  Real PlanckGroupScatteringCoefficientFromNu(const Real rho, const Real temp,
                                              const Real nu) const {
    return ScatteringCoefficientFromNu(rho, temp, nu, Planck);
  }

  PORTABLE_INLINE_FUNCTION
  Real RosselandGroupScatteringCoefficientFromNu(const Real rho,
                                                 const Real temp,
                                                 const Real nu) const {
    return ScatteringCoefficientFromNu(rho, temp, nu, Rosseland);
  }

  PORTABLE_INLINE_FUNCTION
  Real ScatteringCoefficientFromNu(const Real rho, const Real temp,
                                   const Real nu,
                                   const int gmode = Rosseland) const {
    return PortsOfCall::visit(
        [=](const auto &s_opac) {
          return s_opac.ScatteringCoefficientFromNu(rho, temp, nu, gmode);
        },
        s_opac_);
  }

  inline void Finalize() noexcept {
    return PortsOfCall::visit([](auto &s_opac) { return s_opac.Finalize(); },
                              s_opac_);
  }
};

} // namespace impl
} // namespace photons
} // namespace singularity

#endif // SINGULARITY_OPAC_PHOTONS_MEAN_PHOTON_S_VARIANT_
