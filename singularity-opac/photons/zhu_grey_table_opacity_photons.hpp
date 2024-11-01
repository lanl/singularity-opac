// ======================================================================
// © 2024. Triad National Security, LLC. All rights reserved.  This
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

#ifndef SINGULARITY_OPAC_PHOTONS_ZHU_GREY_TABLE_OPACITY_PHOTONS_
#define SINGULARITY_OPAC_PHOTONS_ZHU_GREY_TABLE_OPACITY_PHOTONS_

#include <fstream>
#include <filesystem>

#include <singularity-opac/base/opac_error.hpp>
#include <spiner/databox.hpp>

namespace singularity {
namespace photons {

#define EPS (10.0 * std::numeric_limits<Real>::min())

// Tables have grey opacity dust, molec
// See "Global 3D Radiation Hydrodynamic Simulations of Proto-Jupiter’s
// Convective EnvelopeConvective Envelope" by Zhaohuan Zhu et al (2021)

// TODO: Should the Zhu repository be loaded as a submodule?
// NOTE: Tables are assumed to be log-log in rho-T

// DataBox use is based on, e.g., photon/mean_s_opacity_photons.hpp
template <typename pc = PhysicalConstantsCGS>
class ZhuTableOpacity {
 public:
  using PC = pc;
  using DataBox = Spiner::DataBox<Real>;

  ZhuTableOpacity() = default;

  // construct Planck/Rosseland DataBox from Zhu ascii or spiner hdf5 file
  ZhuTableOpacity(const std::string filename, const bool use_planck_absorb = false)
    : opac_type_(use_planck_absorb) {

    // get number of density and temperature points
    std::ifstream ff(filename.c_str());
    const bool fexists = ff.good();

    if (fexists) {

      std::filesystem::path filePath(filename);
      std::string extension = filePath.extension().string();

      if (extension == ".txt") {
        // assume this is one of the original Zhu et al (2021) ASCII files
        loadZhuASCII(ff);
#ifdef SPINER_USE_HDF
      } else if (extension == ".hdf5" || extension == ".h5") {
        hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        kappa_.loadHDF(file, "Zhu Table");
#endif
      } else {
        OPAC_ERROR("photons::ZhuTableOpacity: unrecognized file extension");
      }

    } else {
      OPAC_ERROR("photons::ZhuTableOpacity: file does not exist");
    }
  }

#ifdef SPINER_USE_HDF
  void Save(const std::string &filename) const {
    herr_t status = H5_SUCCESS;
    hid_t file =
        H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    status += kappa_.saveHDF(file, "Zhu Table");
    status += H5Fclose(file);

    if (status != H5_SUCCESS) {
      OPAC_ERROR("photons::ZhuTableOpacity: HDF5 error\n");
    }
  }
#endif

  ZhuTableOpacity GetOnDevice() { return *this; }
  inline void Finalize() noexcept {}

  PORTABLE_INLINE_FUNCTION
  Real AbsorptionCoefficient(const Real rho, const Real temp, const Real nu,
                             Real *lambda = nullptr) const {
    return dist_.AbsorptionCoefficientFromKirkhoff(*this, rho, temp, nu,
                                                   lambda);
  }

  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void
  AbsorptionCoefficient(const Real rho, const Real temp,
                        FrequencyIndexer &nu_bins, DataIndexer &coeffs,
                        const int nbins, Real *lambda = nullptr) const {
    for (int i = 0; i < nbins; ++i) {
      coeffs[i] = AbsorptionCoefficient(rho, temp, nu_bins[i], lambda);
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real AngleAveragedAbsorptionCoefficient(const Real rho, const Real temp,
                                          const Real nu,
                                          Real *lambda = nullptr) const {
    return dist_.AngleAveragedAbsorptionCoefficientFromKirkhoff(
        *this, rho, temp, nu, lambda);
  }

  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void AngleAveragedAbsorptionCoefficient(
      const Real rho, const Real temp, FrequencyIndexer &nu_bins,
      DataIndexer &coeffs, const int nbins, Real *lambda = nullptr) const {
    for (int i = 0; i < nbins; ++i) {
      coeffs[i] =
          AngleAveragedAbsorptionCoefficient(rho, temp, nu_bins[i], lambda);
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real EmissivityPerNuOmega(const Real rho, const Real temp, const Real nu,
                            Real *lambda = nullptr) const {
    const Real lRho = toLog_(rho);
    const Real lT = toLog_(temp);
    const Real Bnu = dist_.ThermalDistributionOfTNu(temp, nu, lambda);
    return rho * kappa_.interpToReal(lRho, lT, opac_type_) * Bnu;
  }

  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void
  EmissivityPerNuOmega(const Real rho, const Real temp,
                       FrequencyIndexer &nu_bins, DataIndexer &coeffs,
                       const int nbins, Real *lambda = nullptr) const {
    for (int i = 0; i < nbins; ++i) {
      coeffs[i] = EmissivityPerNuOmega(rho, temp, nu_bins[i], lambda);
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real EmissivityPerNu(const Real rho, const Real temp, const Real nu,
                       Real *lambda = nullptr) const {
    return 4 * M_PI * EmissivityPerNuOmega(rho, temp, nu, lambda);
  }

  template <typename FrequencyIndexer, typename DataIndexer>
  PORTABLE_INLINE_FUNCTION void
  EmissivityPerNu(const Real rho, const Real temp, FrequencyIndexer &nu_bins,
                  DataIndexer &coeffs, const int nbins,
                  Real *lambda = nullptr) const {
    for (int i = 0; i < nbins; ++i) {
      coeffs[i] = EmissivityPerNu(rho, temp, nu_bins[i], lambda);
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real Emissivity(const Real rho, const Real temp,
                  Real *lambda = nullptr) const {
    const Real lRho = toLog_(rho);
    const Real lT = toLog_(temp);
    const Real B = dist_.ThermalDistributionOfT(temp, lambda);
    return rho * kappa_.interpToReal(lRho, lT, opac_type_) * B;
  }

  PORTABLE_INLINE_FUNCTION
  Real NumberEmissivity(const Real rho, const Real temp,
                        Real *lambda = nullptr) const {
    const Real lRho = toLog_(rho);
    const Real lT = toLog_(temp);
    return kappa_.interpToReal(lRho, lT, opac_type_) * dist_.ThermalNumberDistributionOfT(temp, lambda);
  }

  PORTABLE_INLINE_FUNCTION
  Real ThermalDistributionOfTNu(const Real temp, const Real nu,
                                Real *lambda = nullptr) const {
    return dist_.ThermalDistributionOfTNu(temp, nu, lambda);
  }

  PORTABLE_INLINE_FUNCTION
  Real DThermalDistributionOfTNuDT(const Real temp, const Real nu,
                                   Real *lambda = nullptr) const {
    return dist_.DThermalDistributionOfTNuDT(temp, nu, lambda);
  }

  PORTABLE_INLINE_FUNCTION
  Real ThermalDistributionOfT(const Real temp, Real *lambda = nullptr) const {
    return dist_.ThermalDistributionOfT(temp, lambda);
  }

  PORTABLE_INLINE_FUNCTION Real
  ThermalNumberDistributionOfT(const Real temp, Real *lambda = nullptr) const {
    return dist_.ThermalNumberDistributionOfT(temp, lambda);
  }

  PORTABLE_INLINE_FUNCTION
  Real EnergyDensityFromTemperature(const Real temp,
                                    Real *lambda = nullptr) const {
    return dist_.EnergyDensityFromTemperature(temp, lambda);
  }

  PORTABLE_INLINE_FUNCTION
  Real TemperatureFromEnergyDensity(const Real er,
                                    Real *lambda = nullptr) const {
    return dist_.TemperatureFromEnergyDensity(er, lambda);
  }

  PORTABLE_INLINE_FUNCTION
  Real NumberDensityFromTemperature(const Real temp,
                                    Real *lambda = nullptr) const {
    return dist_.NumberDensityFromTemperature(temp, lambda);
  }

 private:
  // stolen from mean_s_opacity_photons.hpp
  PORTABLE_INLINE_FUNCTION Real toLog_(const Real x) const {
    return std::log10(std::abs(x) + EPS);
  }
  PORTABLE_INLINE_FUNCTION Real fromLog_(const Real lx) const {
    return std::pow(10., lx);
  }

  // if .txt file, assume it is the original Zhu dust opacity file
  void loadZhuASCII(std::ifstream &ff) {

    int NRho = -1;
    int NT = -1;

    // line read from file
    std::string fline;

    // read 1-line header to get sizes
    std::getline(ff, fline);

    // tokenize fline
    char* cfline = const_cast<char*>(fline.c_str());
    char* fl_tok = std::strtok(cfline, " ");

    // move to next token to get number of density points
    fl_tok = std::strtok(nullptr, " ");
    NRho = std::stoi(fl_tok);

    // move to next token to get number of density points
    fl_tok = std::strtok(nullptr, " ");
    fl_tok = std::strtok(nullptr, " ");
    NT = std::stoi(fl_tok);

    // reseize the Planck and Rosseland databoxes (number of types of opac=2)
    kappa_.resize(NRho, NT, 2);

    // set rho-T rankges (Zhu tables are uniform in log-log rho-T space)
    const Real lTMin = toLog_(1.0);
    const Real lTMax = toLog_(7943282.347242886);
    const Real lRhoMin = toLog_(1.0e-14);
    const Real lRhoMax = toLog_(0.7943282347241912);
    kappa_.setRange(1, lTMin, lTMax, NT);
    kappa_.setRange(2, lRhoMin, lRhoMax, NRho);

    // fill tables
    for (int iRho = 0; iRho < NRho; ++iRho) {
      const Real lRho_i = kappa_.range(2).x(iRho);
      for (int iT = 0; iT < NT; ++iT) {

        // get new file-line
        std::getline(ff, fline);
        cfline = const_cast<char*>(fline.c_str());
        fl_tok = std::strtok(cfline, " ");

        // check for consistent density [g/cm^3] on table row
        const Real Rho = std::stod(fl_tok);
        if (std::abs(Rho - fromLog_(lRho_i)) > 1e-6 * std::abs(Rho)) {
          OPAC_ERROR("photons::ZhuTableOpacity: invalid rho");
        }

        // check for consistent temperature [K] on table row
        const Real lT_i = kappa_.range(1).x(iT);
        fl_tok = std::strtok(nullptr, " ");
        const Real T = std::stod(fl_tok);
        if (std::abs(T - fromLog_(lT_i)) > 1e-6 * std::abs(T)) {
          OPAC_ERROR("photons::ZhuTableOpacity: invalid T");
        }

        // populate Rosseland opacity [cm^2/g]
        fl_tok = std::strtok(nullptr, " ");
        kappa_(iRho, iT, 0) = std::stod(fl_tok);

        // populate Planck opacity [cm^2/g]
        fl_tok = std::strtok(nullptr, " ");
        kappa_(iRho, iT, 1) = std::stod(fl_tok);
      }
    }
  }

  // WARNING: which only one of the two possibilities can be used for now
  int opac_type_;
  PlanckDistribution<pc> dist_;
  DataBox kappa_;
};

} // namespace photons
} // namespace singularity

#endif // SINGULARITY_OPAC_PHOTONS_ZHU_GREY_TABLE_OPACITY_PHOTONS_
