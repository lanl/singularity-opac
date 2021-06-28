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

#ifndef SINGULARITY_OPAC_NEUTRINOS_SPINER_OPAC_NEUTRINOS_HPP_
#define SINGULARITY_OPAC_NEUTRINOS_SPINER_OPAC_NEUTRINOS_HPP_

#include <cassert>
#include <cmath>
#include <cstdio>

#include <ports-of-call/portability.hpp>
#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>

#include <singularity-opac/base/opac_error.hpp>
#include <singularity-opac/base/radiation_types.hpp>
#include <singularity-opac/constants/constants.hpp>

// JMM: Doing everything in log-log, because everything should be
// positive and should roughly follow a power law actually.

namespace singularity {
namespace neutrinos {

namespace impl {
struct TableBounds {
  Real min, max;
  int N;
};
enum class DataStatus { Deallocated, OnDevice, OnHost };
} // namespace impl

class SpinerOpacity {
 public:
  static constexpr Real EPS = 10.0 * std::numeric_limits<Real>::epsilon();

  SpinerOpacity() = default;
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return 0; }
  PORTABLE_INLINE_FUNCTION
  void PrintParams() const noexcept {
    printf("Spiner opacity\n"); // TODO(JMM): Params
  }

  SpinerOpacity GetOnDevice() {
    SpinerOpacity other;
    other.lalphanu_ = Spiner::getOnDeviceDataBox(lalphanu_);
    other.ljnu_ = Spiner::getOnDeviceDataBox(ljnu_);
    other.lJ_ = Spiner::getOnDeviceDataBox(lJ_);
    other.lJYe_ = Spiner::getOnDeviceDataBox(lJYe_);
    other.lRhoBounds_ = lRhoBounds_;
    other.lTBounds_ = lTBounds_;
    other.YeBounds_ = YeBounds_;
    other.lnuBounds_ = lnuBounds_;
    other.memoryStatus_ = impl::DataStatus::OnDevice;
    return other;
  }

  void Finalize() {
    lalphanu_.finalize();
    ljnu_.finalize();
    lJ_.finalize();
    lJYe_.finalize();
  }

 private:
  // TODO(JMM): Offsets probably not necessary
  PORTABLE_INLINE_FUNCTION Real toLog_(const Real x, const Real offset) const {
    return BDMath::log10(std::abs(std::max(x, -offset) + offset) + EPS);
  }
  PORTABLE_INLINE_FUNCTION Real toLog_(const Real x) const {
    return toLog_(x, 0);
  }
  PORTABLE_INLINE_FUNCTION Real fromLog_(const Real lx,
                                         const Real offset) const {
    return std::pow(10., lx) - offset;
  }
  PORTABLE_INLINE_FUNCTION Real fromLog_(const Real lx) const {
    return fromLog_(lx, 0);
  }
  const char *filename_;
  impl::DataStatus memoryStatus_;
  impl::TableBounds lRhoBounds_, lTBounds_, YeBounds_, lnuBounds_;
  // TODO(JMM): Integrating J and JYe seems wise.
  // We can add more things here as needed.
  Spiner::DataBox lalphanu_, ljnu_, lJ_, lJYe_;
};

} // namespace neutrinos
} // namespace singularity

#endif //  SINGULARITY_OPAC_NEUTRINOS_SPINER_OPAC_NEUTRINOS_HPP_
