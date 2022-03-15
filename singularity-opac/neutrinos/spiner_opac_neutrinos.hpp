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
#include <string>

#include <fast-math/logs.hpp>
#include <ports-of-call/portability.hpp>
#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>
#include <spiner/spiner_types.hpp>

#include <singularity-opac/base/opac_error.hpp>
#include <singularity-opac/base/radiation_types.hpp>
#include <singularity-opac/base/sp5.hpp>
#include <singularity-opac/constants/constants.hpp>

#ifdef SPINER_USE_HDF
#include "hdf5.h"
#include "hdf5_hl.h"
#endif

// JMM: Doing everything in log-log, because everything should be
// positive and should roughly follow a power law actually.

// TODO(JMM): The dynamic range in the mantissa is too large, so we
// can't use log10 for floats and take advantage of that performance
// enhancement, unfortunately.
// We should experiment with logs further to see if there's some games
// we can play, but for now, I'm switching everything in this library
// to log10.

namespace singularity {
namespace neutrinos {

namespace impl {
enum class DataStatus { Deallocated, OnDevice, OnHost };
} // namespace impl

// TODO(JMM): Bottom of the table and top of the table handled by
// DataBox. Bottom of the table is a floor. Top of the table is
// power law extrapolation.
// TODO(JMM): Should lJ be stored on disk or computed at start up?
template <typename ThermalDistribution>
class SpinerOpacity {
 public:
  static constexpr Real EPS = 10.0 * std::numeric_limits<Real>::epsilon();
  using pc = PhysicalConstants<CGS>;
  static constexpr Real Hz2MeV = pc::h / (1e6 * pc::eV);
  static constexpr Real MeV2Hz = 1 / Hz2MeV;

  SpinerOpacity() = default;

  // Testing constructor that fills the tables with gray opacities
  template <typename Opacity>
  SpinerOpacity(Opacity &opac, Real lRhoMin, Real lRhoMax, int NRho, Real lTMin,
                Real lTMax, int NT, Real YeMin, Real YeMax, int NYe, Real leMin,
                Real leMax, int Ne)
      : filename_("none"), memoryStatus_(impl::DataStatus::OnHost) {
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
    lJYe_.copyMetadata(lJ_);

    // Fill tables
    for (int iRho = 0; iRho < NRho; ++iRho) {
      Real lRho = lalphanu_.range(4).x(iRho);
      Real rho = fromLog_(lRho);
      for (int iT = 0; iT < NT; ++iT) {
        Real lT = lalphanu_.range(3).x(iT);
        Real T = fromLog_(lT);
        for (int iYe = 0; iYe < NYe; ++iYe) {
          Real Ye = lalphanu_.range(2).x(iYe);
          for (int idx = 0; idx < NEUTRINO_NTYPES; ++idx) {
            RadiationType type = Idx2RadType(idx);
            Real J = std::max(opac.Emissivity(rho, T, Ye, type), 0.0);
            Real lJ = toLog_(J);
            lJ_(iRho, iT, iYe, idx) = lJ;
            Real JYe = std::max(opac.NumberEmissivity(rho, T, Ye, type), 0.0);
            lJYe_(iRho, iT, iYe, idx) = toLog_(JYe);
            for (int ie = 0; ie < Ne; ++ie) {
              Real lE = lalphanu_.range(0).x(ie);
              Real E = fromLog_(lE);
              Real nu = MeV2Hz * E;
              Real alpha = std::max(
                  opac.AbsorptionCoefficient(rho, T, Ye, type, nu), 0.0);
              lalphanu_(iRho, iT, iYe, idx, ie) = toLog_(alpha);
              Real j = std::max(opac.EmissivityPerNuOmega(rho, T, Ye, type, nu),
                                0.0);
              ljnu_(iRho, iT, iYe, idx, ie) = toLog_(j);
            }
          }
        }
      }
    }
  }

  // DataBox constructor. Note that this constructor *shallow* copies
  // the databoxes, so they must be managed externally.
  SpinerOpacity(const Spiner::DataBox &lalphanu, const Spiner::DataBox ljnu,
                const Spiner::DataBox lJ, const Spiner::DataBox lJYe)
      : memoryStatus_(impl::DataStatus::OnHost), lalphanu_(lalphanu),
        ljnu_(ljnu), lJ_(lJ), lJYe_(lJYe) {}

#ifdef SPINER_USE_HDF
  SpinerOpacity(const std::string &filename)
      : filename_(filename.c_str()), memoryStatus_(impl::DataStatus::OnHost) {
    herr_t status = H5_SUCCESS;
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    status += lalphanu_.loadHDF(file, SP5::Opac::AbsorptionCoefficient);
    status += ljnu_.loadHDF(file, SP5::Opac::EmissivityPerNu);
    status += lJ_.loadHDF(file, SP5::Opac::TotalEmissivity);
    status += lJYe_.loadHDF(file, SP5::Opac::NumberEmissivity);
    status += H5Fclose(file);

    if (status != H5_SUCCESS) {
      OPAC_ERROR("neutrinos::SpinerOpacity: HDF5 error\n");
    }
  }

  void Save(const std::string &filename) const {
    herr_t status = H5_SUCCESS;
    hid_t file =
        H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    status += lalphanu_.saveHDF(file, SP5::Opac::AbsorptionCoefficient);
    status += ljnu_.saveHDF(file, SP5::Opac::EmissivityPerNu);
    status += lJ_.saveHDF(file, SP5::Opac::TotalEmissivity);
    status += lJYe_.saveHDF(file, SP5::Opac::NumberEmissivity);
    status += H5Fclose(file);

    if (status != H5_SUCCESS) {
      OPAC_ERROR("neutrinos::SpinerOpacity: HDF5 error\n");
    }
  }
#endif

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
  Real AbsorptionCoefficient(const Real rho, const Real temp, const Real Ye,
                             const RadiationType type, const Real nu,
                             Real *lambda = nullptr) const {
    int idx;
    Real lRho, lT;
    toLogs_(rho, temp, type, lRho, lT, idx);
    const Real le = toLog_(Hz2MeV * nu);
    const Real lalpha = lalphanu_.interpToReal(lRho, lT, Ye, idx, le);
    const Real alpha = fromLog_(lalpha);
    return alpha;
  }

  // TODO(JMM): Should we provide a raw copy operator instead of or
  // addition to interpolation?
  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void
  AbsorptionCoefficient(const Real rho, const Real temp, const Real Ye,
                        const RadiationType type, FrequencyIndexer &nu_bins,
                        DataIndexer &coeffs, const int nbins,
                        Real *lambda = nullptr) const {
    int idx;
    Real lRho, lT;
    toLogs_(rho, temp, type, lRho, lT, idx);
    for (int i = 0; i < nbins; ++i) {
      const Real le = toLog_(Hz2MeV * nu_bins[i]);
      coeffs[i] = fromLog_(lalphanu_.interpToReal(lRho, lT, Ye, idx, le));
    }
  }

  // Angle-averaged absorption coefficient assumed to be the same as absorption
  // coefficient
  PORTABLE_INLINE_FUNCTION
  Real AngleAveragedAbsorptionCoefficient(const Real rho, const Real temp,
                                          const Real Ye,
                                          const RadiationType type,
                                          const Real nu,
                                          Real *lambda = nullptr) const {
    int idx;
    Real lRho, lT;
    toLogs_(rho, temp, type, lRho, lT, idx);
    const Real le = toLog_(Hz2MeV * nu);
    const Real lAlpha = lalphanu_.interpToReal(lRho, lT, Ye, idx, le);
    const Real Alpha = fromLog_(lAlpha);
    return Alpha;
  }

  // TODO(JMM): Should we provide a raw copy operator instead of or
  // addition to interpolation?
  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void AngleAveragedAbsorptionCoefficient(
      const Real rho, const Real temp, const Real Ye, const RadiationType type,
      FrequencyIndexer &nu_bins, DataIndexer &coeffs, const int nbins,
      Real *lambda = nullptr) const {
    int idx;
    Real lRho, lT;
    toLogs_(rho, temp, type, lRho, lT, idx);
    for (int i = 0; i < nbins; ++i) {
      const Real le = toLog_(Hz2MeV * nu_bins[i]);
      coeffs[i] = fromLog_(lalphanu_.interpToReal(lRho, lT, Ye, idx, le));
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real EmissivityPerNuOmega(const Real rho, const Real temp, const Real Ye,
                            const RadiationType type, const Real nu,
                            Real *lambda = nullptr) const {
    int idx;
    Real lRho, lT;
    toLogs_(rho, temp, type, lRho, lT, idx);
    const Real le = toLog_(Hz2MeV * nu);
    const Real lj = ljnu_.interpToReal(lRho, lT, Ye, idx, le);
    return fromLog_(lj);
  }

  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void
  EmissivityPerNuOmega(const Real rho, const Real temp, const Real Ye,
                       const RadiationType type, FrequencyIndexer &nu_bins,
                       DataIndexer &coeffs, const int nbins,
                       Real *lambda = nullptr) const {
    int idx;
    Real lRho, lT;
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
                  const RadiationType type, FrequencyIndexer &nu_bins,
                  DataIndexer &coeffs, const int nbins,
                  Real *lambda = nullptr) const {
    int idx;
    Real lRho, lT;
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
    int idx;
    Real lRho, lT;
    toLogs_(rho, temp, type, lRho, lT, idx);
    const Real lJ = lJ_.interpToReal(lRho, lT, Ye, idx);
    const Real J = fromLog_(lJ);
    return J;
  }

  PORTABLE_INLINE_FUNCTION
  Real NumberEmissivity(const Real rho, const Real temp, Real Ye,
                        RadiationType type, Real *lambda = nullptr) const {
    int idx;
    Real lRho, lT;
    toLogs_(rho, temp, type, lRho, lT, idx);
    return fromLog_(lJYe_.interpToReal(lRho, lT, Ye, idx));
  }

  PORTABLE_INLINE_FUNCTION
  Real ThermalDistributionOfTNu(const Real temp, const RadiationType type,
                                const Real nu, Real *lambda = nullptr) const {
    return dist_.ThermalDistributionOfTNu(temp, type, nu, lambda);
  }

  PORTABLE_INLINE_FUNCTION
  Real ThermalDistributionOfT(const Real temp, const RadiationType type,
                              Real *lambda = nullptr) const {
    return dist_.ThermalDistributionOfT(temp, type, lambda);
  }

  PORTABLE_INLINE_FUNCTION Real ThermalNumberDistributionOfT(
      const Real temp, const RadiationType type, Real *lambda = nullptr) const {
    return dist_.ThermalNumberDistributionOfT(temp, type, lambda);
  }

 private:
  // TODO(JMM): Offsets probably not necessary
  PORTABLE_INLINE_FUNCTION Real toLog_(const Real x, const Real offset) const {
    return std::log10(std::abs(std::max(x, -offset) + offset) + EPS);
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
                                        Real &lT, int &idx) const {
    lRho = toLog_(rho);
    lT = toLog_(temp);
    idx = RadType2Idx(type);
  }
  const char *filename_;
  impl::DataStatus memoryStatus_ = impl::DataStatus::Deallocated;
  // TODO(JMM): Integrating J and JYe seems wise.
  // We can add more things here as needed.
  Spiner::DataBox lalphanu_, ljnu_, lJ_, lJYe_;
  // TODO(JMM): Should we add table bounds? Given they're recorded in
  // each spiner table, I lean towards no, but could be convinced
  // otherwise if we need to do extrapolation, etc.
  ThermalDistribution dist_;
};

} // namespace neutrinos
} // namespace singularity

#endif //  SINGULARITY_OPAC_NEUTRINOS_SPINER_OPAC_NEUTRINOS_HPP_
