// ======================================================================
// © 2021-2024. Triad National Security, LLC. All rights reserved.  This
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
#include <pair>
#include <string>
#include <vector>

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
template <typename ThermalDistribution, typename pc = PhysicalConstantsCGS>
class SpinerOpacity {
 public:
  using PC = pc;
  using DataBox = Spiner::DataBox<Real>;

  static constexpr Real MEV = 1e6 * pc::eV;
  static constexpr Real EPS = 10.0 * std::numeric_limits<Real>::min();
  static constexpr Real Hz2MeV = pc::h / (1e6 * pc::eV);
  static constexpr Real MeV2Hz = 1 / Hz2MeV;
  static constexpr Real MeV2GK = 11.604525006;
  static constexpr Real GK2MeV = 1. / MeV2GK;
  static constexpr Real MeV2K = 1.e9 * MeV2GK;
  static constexpr Real K2MeV = 1. / MeV2K;

  enum class LoadSource { SP5, NuLib };

  SpinerOpacity() = default;

  // Testing constructor that fills the tables with gray opacities
  template <typename Opacity>
  SpinerOpacity(Opacity &opac, Real lRhoMin, Real lRhoMax, int NRho, Real lTMin,
                Real lTMax, int NT, Real YeMin, Real YeMax, int NYe, Real leMin,
                Real leMax, int Ne)
      : filename_("none"), memoryStatus_(impl::DataStatus::OnHost) {
    lTMin += std::log10(K2MeV);
    lTMax += std::log10(K2MeV);
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
            Real J = std::max(opac.Emissivity(rho, T * MeV2K, Ye, type), 0.0);
            Real lJ = toLog_(J);
            lJ_(iRho, iT, iYe, idx) = lJ;
            Real JYe =
                std::max(opac.NumberEmissivity(rho, T * MeV2K, Ye, type), 0.0);
            lJYe_(iRho, iT, iYe, idx) = toLog_(JYe);
            for (int ie = 0; ie < Ne; ++ie) {
              Real lE = lalphanu_.range(0).x(ie);
              Real E = fromLog_(lE);
              Real nu = MeV2Hz * E;
              Real alpha = std::max(
                  opac.AbsorptionCoefficient(rho, T, Ye, type, nu), 0.0);
              lalphanu_(iRho, iT, iYe, idx, ie) = toLog_(alpha);
              Real j = std::max(
                  opac.EmissivityPerNuOmega(rho, T * MeV2K, Ye, type, nu), 0.0);
              ljnu_(iRho, iT, iYe, idx, ie) = toLog_(j);
            }
          }
        }
      }
    }
  }

  // DataBox constructor. Note that this constructor *shallow* copies
  // the databoxes, so they must be managed externally.
  SpinerOpacity(const DataBox &lalphanu, const DataBox ljnu, const DataBox lJ,
                const DataBox lJYe)
      : memoryStatus_(impl::DataStatus::OnHost), lalphanu_(lalphanu),
        ljnu_(ljnu), lJ_(lJ), lJYe_(lJYe) {}

#ifdef SPINER_USE_HDF
  SpinerOpacity(const std::string &filename,
                LoadSource load_from = LoadSource::SP5)
      : filename_(filename.c_str()), memoryStatus_(impl::DataStatus::OnHost) {
    if (load_from == LoadSource::SP5) {
      LoadFromSP5_(filename);
    } else if (load_from == LoadSource::NuLib) {
      LoadFromNuLib_(filename);
    } else {
      OPAC_ERROR("Unknown file source type\n");
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
#endif // SPINER_USE_HDF

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
  Real DThermalDistributionOfTNuDT(const Real temp, const RadiationType type,
                                   const Real nu,
                                   Real *lambda = nullptr) const {
    return dist_.DThermalDistributionOfTNuDT(temp, type, nu, lambda);
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

  PORTABLE_INLINE_FUNCTION
  Real EnergyDensityFromTemperature(const Real temp, const RadiationType type,
                                    Real *lambda = nullptr) const {
    return dist_.EnergyDensityFromTemperature(temp, type, lambda);
  }

  PORTABLE_INLINE_FUNCTION
  Real TemperatureFromEnergyDensity(const Real er, const RadiationType type,
                                    Real *lambda = nullptr) const {
    return dist_.TemperatureFromEnergyDensity(er, type, lambda);
  }

  PORTABLE_INLINE_FUNCTION
  Real NumberDensityFromTemperature(const Real temp, const RadiationType type,
                                    Real *lambda = nullptr) const {
    return dist_.NumberDensityFromTemperature(temp, type, lambda);
  }

 private:
#ifdef SPINER_USE_HDF
  void LoadFromSP5_(const std::string &filename) {
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

  void LoadFromNuLib_(const std::string &filename) {
    herr_t status = H5_SUCCESS;
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    // table size
    const int NR = ReadInt_(file, "nrho");
    const int NT = ReadInt_(file, "ntemp");
    const int NY = ReadInt_(file, "nye");
    const int NE = ReadInt_(file, "number_groups");
    const int NTYPE = ReadInt_(file, "number_species");

    auto [rho_min, rho_max] = ReadBounds_(file, "rho_points", NR);
    auto [T_min, T_max] = ReadBounds_(file, "temp_points", NT);
    auto [Ye_min, Ye_max] = ReadBounds_(file, "ye_points", NY);
    auto [emin, emax] = ReadBounds_(file, "neutrino_energies", NE);

    const Real lRhoMin = std::log10(rho_min); // g/cm^3
    const Real lRhoMax = std::log10(rho_max);
    const Real lTMin = std::log10(T_min); // MeV internally. Converted from K
    const Real lTMax = std::log10(T_max); // when function calls made.
    const Real leMin = std::log10(emin);  // MeV internally. Converted from Hz
    const Real leMax = std::log10(emax);  // when function calls made.

    DataBox jnu_file(NE, NTYPE, NY, NT,
                     NR); // (erg s^{-1} cm^{-3} MeV^{-1} / 4pi)
    status = H5LTread_dataset_double(file, "emissivities", jnu_file.data());
    if (status != H5_SUCCESS) {
      OPAC_ERROR("An HDF5 error ocurred while reading emissivities");
    }
    DataBox kappa_file(NE, NTYPE, NY, NT, NR); // 1/cm
    status =
        H5LTread_dataset_double(file, "absorption_opacity", jnu_file.data());
    if (status != H5_SUCCESS) {
      OPAC_ERROR("An HDF5 error ocurred while reading absorption opacities");
    }
    status = H5Fclose(file);
    if (status != H5_SUCCESS) {
      OPAC_ERROR("An HDF5 error ocurred while reading absorption opacities");
    }

    // Set metadata for lalphanu and ljnu
    lalphanu_.resize(NR, NT, NY, NTYPE, NE); // 1/cm
    lalphanu_.setRange(0, leMin, leMax, NE);
    // index 1 is the species and is not interpolatable
    lalphanu_.setRange(2, YeMin, YeMax, NY);
    lalphanu_.setRange(3, lTMin, lTMax, NT);
    lalphanu_.setRange(4, lRhoMin, lRhoMax, NR);
    ljnu_.copyMetadata(lalphanu_); // (erg s^{-1} cm^{-3} Hz^{-1} / 4pi)

    // set metadata for lJ and lJYe
    lJ_.resize(NR, NT, NY, NTYPE);
    lJ_.setRange(1, YeMin, YeMax, NY);
    lJ_.setRange(2, lTMin, lTMax, NT);
    lJ_.setRange(3, lRhoMin, lRhoMax, NR);
    lJYe_.copyMetadata(lJ_);

    // Fill lalphanu and lJnu
    for (int iR = 0; iR < NR; ++iR) {
      for (int iT = 0; iT < NT; ++iT) {
        for (int iY = 0; iY < NY, ++iY) {
          for (int idx = 0; idx < NTYPE; ++idx) {
            for (int ie = 0; ie < NE; ++ie) {
              Real kappa = kappa_file(ie, idx, iY, iT, iR);
              // log
              lalphanu_(iR, iT, iY, idx, ie) = std::log10(kappa);
              Real j = jnu_file(ie, idx, iY, iT, iR);
              // convert from
              // (erg s^{-1} cm^{-3} MeV^{-1} / 4pi)
              // to
              // (erg s^{-1} cm^{-3} Hz^{-1} / 4pi)
              ljnu_(iR, iT, iY, idx, ie) = std::log10(pc::h * j / MEV);
            }
          }
        }
      }
    }
    ComputeIntegrals(ljnu_, lJ_, lJYe_);
  }

  static int ReadInt_(const hid_t &file_id, const std::string &name) {
    int data;
    herr_t status = H5LTread_dataset_int(file_id, name.c_str(), &data);
    if (status != H5_SUCCESS) {
      OPAC_ERROR("Failed to read dataset!\n");
    }
    return data;
  }

  static auto ReadBounds_(const hid_t &file_id, const std::string &name,
                          int size) {
    std::vector<Real> table(size);
    herr_t status =
        H5LTread_dataset_double(file_id, name.c_str(), table.data());
    if (status != H5_SUCCESS) {
      OPAC_ERROR("An HDF5 error ocurred while reading bounds");
    }
    Real lo = table[0];
    Real hi = table[size - 1];
    return std::make_pair(lo, hi);
  }

  static void ComputeIntegrals(const DataBox &ljnu, DataBox &lJ,
                               DataBox &lJYe) {
    auto rhoGrid = ljnu.range(4);
    auto TGrid = ljnu.range(3);
    auto YeGrid = ljnu.range(2);
    auto leGrid = ljnu.range(0);

    // integrals performed as Riemann sum in log space
    Real dle = leGrid.dx();
    for (int iRho = 0; iRho < rhoGrid.nPoints(); ++iRho) {
      for (int iT = 0; iT < TGrid.nPoints(); ++iT) {
        for (int iYe = 0; iYe < YeGrid.nPoints(); ++iYe) {
          for (int itp = 0; itp < RAD_NUM_TYPES; ++itp) {
            lJ(iRho, iT, iYe, itp) = 0.0;
            lJYe(iRho, iT, iYe, itp) = 0.0;
            for (int ie = 0; ie < leGrid.nPoints(); ++ie) {
              Real le = leGrid.x(ie);
              Real e = std::pow(10., le);
              Real lj = ljnu(iRho, iT, iYe, itp, ie);
              // convert again to
              // (erg s^{-1} cm^{-3} MeV^{-1} / 4pi)
              // since we integrate in those units
              Real j = std::pow(10., lj) * MEV / pc::h;
              Real integrand = 4 * M_PI * j * e;
              lJ(iRho, iT, iYe, itp) += integrand;
              // divide by energy in ergs for lJYe to get number emissivity
              lJYe(iRho, iT, iYe, itp) += integrand / (MEV * e);
            }
            // multiply by log spacing and take the log
            lJ(iRho, iT, iYe, itp) = ToLog(lJ(iRho, iT, iYe, itp) * dle);
            lJYe(iRho, iT, iYe, itp) = ToLog(lJYe(iRho, iT, iYe, itp) * dle);
          }
        }
      }
    }
  }
#endif // SPINER_USE_HDF

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
    lT = toLog_(temp * K2MeV);
    idx = RadType2Idx(type);
  }
  const char *filename_;
  impl::DataStatus memoryStatus_ = impl::DataStatus::Deallocated;
  // TODO(JMM): Integrating J and JYe seems wise.
  // We can add more things here as needed.
  DataBox lalphanu_, ljnu_, lJ_, lJYe_;
  // TODO(JMM): Should we add table bounds? Given they're recorded in
  // each spiner table, I lean towards no, but could be convinced
  // otherwise if we need to do extrapolation, etc.
  ThermalDistribution dist_;
};

} // namespace neutrinos
} // namespace singularity

#endif //  SINGULARITY_OPAC_NEUTRINOS_SPINER_OPAC_NEUTRINOS_HPP_
