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

#ifndef SINGULARITY_OPAC_NEUTRINOS_NON_CGS_S_NEUTRINOS_
#define SINGULARITY_OPAC_NEUTRINOS_NON_CGS_S_NEUTRINOS_

#include <cassert>
#include <cmath>
#include <cstdio>

#include <ports-of-call/portability.hpp>
#include <singularity-opac/base/opac_error.hpp>

namespace singularity {
namespace neutrinos {

template <typename SOpac>
class NonCGSUnitsS {
 public:
  NonCGSUnitsS() = default;
  NonCGSUnitsS(SOpac &&s_opac, const Real time_unit, const Real mass_unit,
              const Real length_unit, const Real temp_unit)
      : s_opac_(std::forward<SOpac>(s_opac)), time_unit_(time_unit),
        mass_unit_(mass_unit), length_unit_(length_unit), temp_unit_(temp_unit),
        rho_unit_(mass_unit_ / (length_unit_ * length_unit_ * length_unit_)),
        freq_unit_(1. / time_unit_) {}

  auto GetOnDevice() {
    return NonCGSUnitsS<SOpac>(s_opac_.GetOnDevice(), time_unit_, mass_unit_,
                             length_unit_, temp_unit_);
  }
  inline void Finalize() noexcept { s_opac_.Finalize(); }

  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return s_opac_.nlambda(); }

  PORTABLE_INLINE_FUNCTION
  Real ScatteringCoefficient(const Real rho, const Real temp, const Real Ye,
                             const RadiationType type, const Real nu,
                             Real *lambda = nullptr) const {
    const Real alpha = s_opac_.ScatteringCoefficient(
        rho_unit_ * rho, temp_unit_ * temp, Ye, type, nu * freq_unit_, lambda);
    // alpha output in units of 1/cm. Want to convert out of CGS.
    // multiplication by length_unit converts length to cm.
    // division converts length from cm to unit system.
    // thus multiplication converts (1/cm) to unit system.
    return alpha * length_unit_;
  }

  PORTABLE_INLINE_FUNCTION
  Real AngleAveragedScatteringCoefficient(const Real rho, const Real temp,
                                          const Real Ye,
                                          const RadiationType type,
                                          const Real nu,
                                          Real *lambda = nullptr) const {
    const Real alpha = s_opac_.AngleAveragedScatteringCoefficient(
        rho * rho_unit_, temp * temp_unit_, Ye, type, nu * freq_unit_, lambda);
    return alpha * length_unit_;
  }

 private:
  SOpac s_opac_;
  Real time_unit_, mass_unit_, length_unit_, temp_unit_;
  Real rho_unit_, freq_unit_;
};

template <typename MeanSOpac>
class MeanNonCGSUnitsS {
 public:
  MeanNonCGSUnitsS() = default;
  MeanNonCGSUnitsS(MeanSOpac &&mean_s_opac, const Real time_unit,
                  const Real mass_unit, const Real length_unit,
                  const Real temp_unit)
      : mean_s_opac_(std::forward<MeanSOpac>(mean_s_opac)), time_unit_(time_unit),
        mass_unit_(mass_unit), length_unit_(length_unit), temp_unit_(temp_unit),
        rho_unit_(mass_unit_ / (length_unit_ * length_unit_ * length_unit_)) {}

  auto GetOnDevice() {
    return MeanNonCGSUnitsS<MeanSOpac>(mean_s_opac_.GetOnDevice(), time_unit_,
                                     mass_unit_, length_unit_, temp_unit_);
  }
  inline void Finalize() noexcept { mean_s_opac_.Finalize(); }

  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return mean_s_opac_.nlambda(); }

  PORTABLE_INLINE_FUNCTION
  Real PlanckMeanScatteringCoefficient(const Real rho, const Real temp,
                                       const Real Ye,
                                       const RadiationType type) const {
    const Real alpha = mean_s_opac_.PlanckMeanScatteringCoefficient(
        rho_unit_ * rho, temp_unit_ * temp, Ye, type);
    // alpha output in units of 1/cm. Want to convert out of CGS.
    // multiplication by length_unit converts length to cm.
    // division converts length from cm to unit system.
    // thus multiplication converts (1/cm) to unit system.
    return alpha * length_unit_;
  }

  PORTABLE_INLINE_FUNCTION
  Real RosselandMeanScatteringCoefficient(const Real rho, const Real temp,
                                          const Real Ye,
                                          const RadiationType type) const {
    const Real alpha = mean_s_opac_.RosselandMeanScatteringCoefficient(
        rho_unit_ * rho, temp_unit_ * temp, Ye, type);
    // alpha output in units of 1/cm. Want to convert out of CGS.
    // multiplication by length_unit converts length to cm.
    // division converts length from cm to unit system.
    // thus multiplication converts (1/cm) to unit system.
    return alpha * length_unit_;
  }

 private:
  MeanSOpac mean_s_opac_;
  Real time_unit_, mass_unit_, length_unit_, temp_unit_;
  Real rho_unit_;
};

} // namespace neutrinos
} // namespace singularity

#endif // SINGULARITY_OPAC_NEUTRINOS_NON_CGS_S_NEUTRINOS_
