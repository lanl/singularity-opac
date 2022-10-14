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

#ifndef SINGULARITY_OPAC_NEUTRINOS_GRAY_S_OPACITY_NEUTRINOS_
#define SINGULARITY_OPAC_NEUTRINOS_GRAY_S_OPACITY_NEUTRINOS_

#include <cassert>
#include <cmath>
#include <cstdio>

#include <ports-of-call/portability.hpp>
#include <singularity-opac/base/opac_error.hpp>

namespace singularity {
namespace neutrinos {

template <typename pc = PhysicalConstantsCGS>
class GraySOpacity {
 public:
  GraySOpacity() = default;
  GraySOpacity(const Real sigma, const Real mmw) : sigma_(sigma), mmw_(mmw) {}

  GraySOpacity GetOnDevice() { return *this; }
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return 0; }
  PORTABLE_INLINE_FUNCTION
  void PrintParams() const noexcept {
    printf("Gray scattering opacity. sigma = %g mmw = %g\n", sigma_, mmw_);
  }
  inline void Finalize() noexcept {}

  PORTABLE_INLINE_FUNCTION
  Real TotalCrossSection(const Real rho, const Real temp, const Real Ye,
                         const RadiationType type, const Real nu, Real *lambda = nullptr) const {
                         return sigma_;
                         }

  PORTABLE_INLINE_FUNCTION
  Real DifferentialCrossSection(const Real rho, const Real temp, const Real Ye,
      const RadiationType type, const Real nu, const Real theta, Real *lambda = nullptr) const {
        return sigma_ / (4. * M_PI);
      }

  PORTABLE_INLINE_FUNCTION
  Real TotalScatteringCoefficient(const Real rho, const Real temp, const Real Ye,
                             const RadiationType type, const Real nu,
                             Real *lambda = nullptr) const {
    return (rho / mmw_) * sigma_;
  }

 private:
  Real sigma_; // Scattering cross section. Units of cm^2
  Real mmw_; // Mean molecular weight. Units of g
};

} // namespace neutrinos
} // namespace singularity

#endif // SINGULARITY_OPAC_NEUTRINOS_GRAY_S_OPACITY_NEUTRINOS_
