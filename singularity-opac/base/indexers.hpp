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

#ifndef SINGULARITY_OPAC_BASE_INDEXERS_
#define SINGULARITY_OPAC_BASE_INDEXERS_

#include <fast-math/logs.hpp>
#include <spiner/databox.hpp>
#include <variant/include/mpark/variant.hpp>

// Indexers are filled by by the the opacity functions as a way to
// fill multiple frequency bins/groups at once. All an indexer *must*
// have is an operator[] method that works on device. However, here we
// show that indexers can have multiple other capabilities, such as
// interpolation via a Spiner DataBox.

namespace singularity {
namespace indexers {

class Linear {
 public:
  Linear() = default;

  PORTABLE_INLINE_FUNCTION
  Linear(Real *data, Real numin, Real numax, int N) : data_(data, N) {
    SetRange_(numin, numax, N);
  }

  Linear(Real numin, Real numax, int N)
      : data_(Spiner::AllocationTarget::Device, N) {
    SetRange_(numin, numax, N);
  }

  PORTABLE_INLINE_FUNCTION Linear(const Spiner::DataBox &data, Real numin,
                                  Real numax, int N)
      : data_(data) {
    SetRange_(numin, numax, N);
  }

  void Finalize() { data_.finalize(); }

  PORTABLE_INLINE_FUNCTION
  Real &operator[](const int i) { return data_(i); }
  PORTABLE_INLINE_FUNCTION
  Real &operator[](const int i) const { return data_(i); }

  PORTABLE_INLINE_FUNCTION
  Real operator()(const Real nu) { return data_.interpToReal(nu); }

 private:
  PORTABLE_INLINE_FUNCTION
  void SetRange_(Real numin, Real numax, int N) {
    data_.setRange(0, numin, numax, N);
  }
  Spiner::DataBox data_;
};

class LogLinear {
 public:
  LogLinear() = default;

  PORTABLE_INLINE_FUNCTION
  LogLinear(Real *data, Real numin, Real numax, int N) : data_(data, N) {
    SetRange_(numin, numax, N);
  }

  LogLinear(Real numin, Real numax, int N)
      : data_(Spiner::AllocationTarget::Device, N) {
    SetRange_(numin, numax, N);
  }

  PORTABLE_INLINE_FUNCTION
  LogLinear(const Spiner::DataBox &data, Real numin, Real numax, int N)
      : data_(data) {
    SetRange_(numin, numax, N);
  }

  void Finalize() { data_.finalize(); }

  PORTABLE_INLINE_FUNCTION
  Real &operator[](const int i) { return data_(i); }
  PORTABLE_INLINE_FUNCTION
  Real &operator[](const int i) const { return data_(i); }

  PORTABLE_INLINE_FUNCTION
  Real operator()(const Real nu) {
    return data_.interpToReal(BDMath::log10(nu));
  }

 private:
  PORTABLE_INLINE_FUNCTION
  void SetRange_(Real numin, Real numax, int N) {
    data_.setRange(0, BDMath::log10(numin), BDMath::log10(numax), N);
  }
  Spiner::DataBox data_;
};

} // namespace indexers
} // namespace singularity

#endif // SINGULARITY_BASE_NEUTRINOS_INDEXERS_
