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

#include <fast-math/logs.hpp>
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
  TableBounds(const Real &min_, const Real &max_, const int &N_)
      : min_(min), max_(max), N_(N);
  Real min, max;
  int N;
};
enum class DataStatus { Deallocated, OnDevice, OnHost };
} // namespace impl

// TODO(JMM): Bottom of the table and top of the table handled by
// DataBox. Bottom of the table is a floor. Top of the table is
// power law extrapolation.
class SpinerOpacity {
 public:
  static constexpr Real EPS = 10.0 * std::numeric_limits<Real>::epsilon();
  constexpr PhysicalConstants<CGS> pc;
  constexpr Real Hz2MeV = pc.h / (1e6 * pc.eV);
  constexpr Real MeV2Hz = 1 / Hz2Mev;

  SpinerOpacity() = default;

  // Testing constructor that fills the tables with gray opacities
  template <typename Opacity>
  SpinerOpacity(Opacity &opac, Real lRhoMin, Real lRhoMax, int NRho, Real lTMin,
                Real lTMax, int NT, Real YeMin, Real YeMax, int NYe, Real leMin,
                Real leMax, int Ne)
      : lRhoBounds_(lRhoMin, lRhoMax, NRho), lTBounds_(lTMin, lTMax, NT),
        YeBounds(YeMin, YeMax, NYe), leBounds_(leMin, leMax, Ne) {
    // Set metadata for lalphanu and ljnu
    lalphanu_.resize(NRho, NT, NYe, NEUTRINO_NTYPES, Ne);
    lalphanu_.setRange(0, leMin, leMax, Ne);
    // index 1 is the species and is not interpolatable
    lalphanu_.setRange(2, YeMin, YeMax, NYe);
    lalphanu_.setRange(3, lTMin, lTMax, NT);
    lalphanu_.setRange(4, lRhoMin, lRhoMax, NRho);
    ljnu_.copyMetadata(lalphanu_);

    // set metadata for lJ and lJYe
    lJ_.resize(NRho, NT, NYe, NEUTRINO_NTYPES);
    lJ_.setRange(1, YeMin, YeMax, NYe);
    lJ_.setRange(2, lTMin, lTMax, NT);
    lJ_.setRange(3, lRhoMin, lRhoMax, NRho);
    lJYe.copyMetadata(lJ_);

    // Fill tables
    for (int iRho = 0; iRho < NRho; ++iRho) {
      Real lRho = lalphanu_.range(4).x(iRho);
      Real rho = fromLog_(lRho);
      for (int iT = 0; iT < NT; ++iT) {
        Real lT = lalphanu_.range(3).x(iT);
        Real T = fromLog_(lT);
        for (int iYe = 0; iYe < NYe; ++iYe) {
          Real Ye = lalphanu_.range(2).x(iYe);
          for (int type = 0; type < NEUTRINO_NTYPES; ++type) {
            lJ_.(iRho, iT, iYe, type) =
                toLog_(opac.Emissivity(rho, T, Ye, type));
            lJYe_.(iRho, iT, iYe, type) =
                toLog_(opac.Emissivity(rho, T, Ye, type));
            for (int ie = 0; ie < Ne; ++ie) {
              Real lE = lalphanu_.range(0).x(ie);
              Real E = fromLog_(lE);
              Real nu = MeV2Hz * E;
              lalphanu_(iRho, iT, iYe, type, ie) =
                  toLog_(opac.AbsorptionCoefficientPerNu(rho, T, Ye, type, nu));
              ljnu_(iRho, iT, iYe, type, ie) =
                  toLog_(opac.EmissivityPerNuOmega(rho, T, Ye, type, nu));
            }
          }
        }
      }
    }
  }

  // TODO(JMM): Constructor from file
  // TODO(JMM): Should lJ be stored on disk or computed at start up?

  PORTABLE_INLINE_FUNCTION int nlambda() const noexcept { return 0; }
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

  PORTABLE_INLINE_FUNCTION
  Real AbsorptionCoefficientPerNu(const Real rho, const Real temp,
                                  const Rela Ye, const RadiationType type,
                                  const Real nu, Real *lambda = nullptr) const {
    Real lRho, lT, idx;
    toLogs(rho, temp, type, lRho, lT, idx);
    const Real le = toLog_(Hz2MeV * nu);
    return fromLog_(lalphanu_.interpToReal(lRho, lT, Ye, idx, le));
  }

  // TODO(JMM): Should we provide a raw copy operator instead of or
  // addition to interpolation?
  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void AbsportionCoefficientPerNu(
      const Real rho, const Real temp, const Real Ye, const RadiationType type,
      const FrequencyIndexer &nu_bins, DataIndexer &coeffs, const int nbins,
      Real *lambda = nullptr) const {
    Real lRho, lT, idx;
    toLogs_(rho, temp, type, lRho, lT, idx);
    for (int i = 0; i < nbins; ++i) {
      const Real le = toLog_(Hz2MeV * nu_bins[i]);
      coeffs[i] = fromLog_(lalphanu_.interpToReal(lRho, lT, Ye, idx, le));
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real EmissivityPerNuOmega(const Real rho, const Real temp, const Real Ye,
                            const RadiationType type, const Real nu,
                            Real *lambda = nullptr) {
    Real lRho, lT, idx;
    toLogs_(rho, temp, type, lRho, lT, idx);
    const Real le = toLog_(Hz2MeV * nu);
    return fromLog_(ljnu_.interpToReal(lRho, lT, le, idx));
  }

  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void
  EmissivityPerNuOmega(const Real rho, const Real temp, const Real Ye,
                       const RadiationType type,
                       const FrequencyIndexer &nu_bins, DataIndexer &coeffs,
                       const int nbins, Real *lambda = nullptr) const {
    Real lRho, lT, idx;
    toLogs_(rho, temp, type, lRho, lT, idx);
    for (int i = 0; i < nbins; ++i) {
      const Real le = toLog_(Hz2MeV * nu_bins[i]);
      coeffs[i] = fromLog_(ljnu_.interpToReal(lRho, lT, Ye, idx, le));
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real EmissivityPerNu(const Real rho, const Real temp, const Real Ye,
                       const RadiationType type, const Real nu,
                       Real *lambda = nullptr) const {
    return 4 * M_PI * EmissivityPerNuOmega(rho, temp, Ye, type, nu, lambda);
  }

  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void
  EmissivityPerNu(const Real rho, const Real temp, const Real Ye,
                  const RadiationType type, const FrequencyIndexer &nu_bins,
                  DataIndexer &coeffs, const int nbins,
                  Real *lambda = nullptr) const {
    Real lRho, lT, idx;
    toLogs_(rho, temp, type, lRho, lT, idx);
    for (int i = 0; i < nbins; ++i) {
      const Real le = toLog_(Hz2MeV * nu_bins[i]);
      coeffs[i] =
          4 * M_PI * fromLog_(ljnu_.interpToReal(lRho, lT, Ye, idx, le));
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real Emissivity(const Real rho, const Real temp, const Real Ye,
                  const RadiationType type, Real *lambda = nullptr) const {
    Real lRho, lT, idx;
    toLogs_(rho, temp, type, lRho, lT, idx);
    return fromLog_(lJ_.interpToReal(lRho, lT, Ye, idx));
  }

  PORTABLE_INLINE_FUNCTION
  Real NumberEmissivity(const Real rho, const Real temp, Real Ye,
                        RadiationType type, Real *lambda = nullptr) const {
    Real lRho, lT, idx;
    toLogs_(rho, temp, type, lRho, lT, idx);
    return fromLog_(lJYe_.interpToReal(lRho, lT, Ye, idx));
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
  PORTABLE_INLINE_FUNCTION void toLogs_(const Real rho, const Real temp,
                                        RadiationType type, Real &lRho,
                                        Real &lT, int &idx) {
    lRho = toLog_(rho);
    lT = toLog_(temp);
    idx = RadType2Idx(type);
  }
  const char *filename_;
  impl::DataStatus memoryStatus_;
  // TODO(JMM): Should the table be tabulated in frequencies or MeV?
  // MeV probably makes more sense, even if the host code uses
  // frequencies.
  impl::TableBounds lRhoBounds_, lTBounds_, YeBounds_, leBounds_;
  // TODO(JMM): Integrating J and JYe seems wise.
  // We can add more things here as needed.
  Spiner::DataBox lalphanu_, ljnu_, lJ_, lJYe;
};

} // namespace neutrinos
} // namespace singularity

#endif //  SINGULARITY_OPAC_NEUTRINOS_SPINER_OPAC_NEUTRINOS_HPP_
