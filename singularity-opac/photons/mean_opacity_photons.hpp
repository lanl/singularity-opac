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

#ifndef SINGULARITY_OPAC_PHOTONS_MEAN_OPACITY_PHOTONS_
#define SINGULARITY_OPAC_PHOTONS_MEAN_OPACITY_PHOTONS_

#include <cmath>
#include <filesystem>
#include <fstream>

#include <ports-of-call/portability.hpp>
#include <singularity-opac/base/radiation_types.hpp>
#include <singularity-opac/base/robust_utils.hpp>
#include <singularity-opac/base/sp5.hpp>
#include <singularity-opac/constants/constants.hpp>
#include <spiner/databox.hpp>

#include <singularity-opac/photons/mean_photon_variant.hpp>
#include <singularity-opac/photons/non_cgs_photons.hpp>

namespace singularity {
namespace photons {
namespace impl {

#define EPS (10.0 * std::numeric_limits<Real>::min())

// TODO(BRR) Note: It is assumed that lambda is constant for all densities and
// temperatures

class MeanOpacity {

 public:
  MeanOpacity() = default;
  template <typename Opacity>
  MeanOpacity(const Opacity &opac, const Real lRhoMin, const Real lRhoMax,
              const int NRho, const Real lTMin, const Real lTMax, const int NT,
              Real *lambda = nullptr) {
    MeanOpacityImpl_<Opacity, true>(opac, lRhoMin, lRhoMax, NRho, lTMin, lTMax,
                                    NT, -1., -1., 100, lambda);
  }

  template <typename Opacity>
  MeanOpacity(const Opacity &opac, const Real lRhoMin, const Real lRhoMax,
              const int NRho, const Real lTMin, const Real lTMax, const int NT,
              Real lNuMin, Real lNuMax, const int NNu, Real *lambda = nullptr) {
    MeanOpacityImpl_<Opacity, false>(opac, lRhoMin, lRhoMax, NRho, lTMin, lTMax,
                                     NT, lNuMin, lNuMax, NNu, lambda);
  }

  // construct Planck/Rosseland DataBox from Zhu-formatted ascii file
  MeanOpacity(const std::string &filename)  : filename_(filename.c_str()) {

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
      } else if (extension == ".hdf5" || extension == ".h5" || extension == ".sp5") {
        herr_t status = H5_SUCCESS;
        hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        status += lkappaPlanck_.loadHDF(file, SP5::MeanOpac::PlanckMeanOpacity);
        status +=
          lkappaRosseland_.loadHDF(file, SP5::MeanOpac::RosselandMeanOpacity);
        status += H5Fclose(file);

        if (status != H5_SUCCESS) {
          OPAC_ERROR("photons::MeanOpacity: HDF5 error\n");
        }
#endif
      } else {
        OPAC_ERROR("photons::MeanOpacity: unrecognized file extension");
      }

    } else {
      OPAC_ERROR("photons::MeanOpacity: file does not exist");
    }
  }

#ifdef SPINER_USE_HDF
  void Save(const std::string &filename) const {
    herr_t status = H5_SUCCESS;
    hid_t file =
        H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    status += lkappaPlanck_.saveHDF(file, SP5::MeanOpac::PlanckMeanOpacity);
    status +=
        lkappaRosseland_.saveHDF(file, SP5::MeanOpac::RosselandMeanOpacity);
    status += H5Fclose(file);

    if (status != H5_SUCCESS) {
      OPAC_ERROR("photons::MeanOpacity: HDF5 error\n");
    }
  }
#endif

  PORTABLE_INLINE_FUNCTION void PrintParams() const {
    printf("Mean opacity\n");
  }

  MeanOpacity GetOnDevice() {
    MeanOpacity other;
    other.lkappaPlanck_ = Spiner::getOnDeviceDataBox(lkappaPlanck_);
    other.lkappaRosseland_ = Spiner::getOnDeviceDataBox(lkappaRosseland_);
    return other;
  }

  void Finalize() {
    lkappaPlanck_.finalize();
    lkappaRosseland_.finalize();
  }

  PORTABLE_INLINE_FUNCTION
  Real PlanckMeanAbsorptionCoefficient(const Real rho, const Real temp) const {
    Real lRho = toLog_(rho);
    Real lT = toLog_(temp);
    return rho * fromLog_(lkappaPlanck_.interpToReal(lRho, lT));
  }

  PORTABLE_INLINE_FUNCTION
  Real RosselandMeanAbsorptionCoefficient(const Real rho,
                                          const Real temp) const {
    Real lRho = toLog_(rho);
    Real lT = toLog_(temp);
    return rho * fromLog_(lkappaRosseland_.interpToReal(lRho, lT));
  }

 private:
  template <typename Opacity, bool AUTOFREQ>
  void MeanOpacityImpl_(const Opacity &opac, const Real lRhoMin,
                        const Real lRhoMax, const int NRho, const Real lTMin,
                        const Real lTMax, const int NT, Real lNuMin,
                        Real lNuMax, const int NNu, Real *lambda = nullptr) {
    using PC = typename Opacity::PC;

    lkappaPlanck_.resize(NRho, NT);
    lkappaPlanck_.setRange(0, lTMin, lTMax, NT);
    lkappaPlanck_.setRange(1, lRhoMin, lRhoMax, NRho);
    lkappaRosseland_.copyMetadata(lkappaPlanck_);

    // Fill tables
    for (int iRho = 0; iRho < NRho; ++iRho) {
      Real lRho = lkappaPlanck_.range(1).x(iRho);
      Real rho = fromLog_(lRho);
      for (int iT = 0; iT < NT; ++iT) {
        Real lT = lkappaPlanck_.range(0).x(iT);
        Real T = fromLog_(lT);
        Real kappaPlanckNum = 0.;
        Real kappaPlanckDenom = 0.;
        Real kappaRosselandNum = 0.;
        Real kappaRosselandDenom = 0.;
        if (AUTOFREQ) {
          lNuMin = toLog_(1.e-3 * PC::kb * fromLog_(lTMin) / PC::h);
          lNuMax = toLog_(1.e3 * PC::kb * fromLog_(lTMax) / PC::h);
        }
        const Real dlnu = (lNuMax - lNuMin) / (NNu - 1);
        // Integrate over frequency
        for (int inu = 0; inu < NNu; ++inu) {
          const Real weight =
              (inu == 0 || inu == NNu - 1) ? 0.5 : 1.; // Trapezoidal rule
          const Real lnu = lNuMin + inu * dlnu;
          const Real nu = fromLog_(lnu);
          const Real alpha = opac.AbsorptionCoefficient(rho, T, nu, lambda);
          const Real B = opac.ThermalDistributionOfTNu(T, nu);
          const Real dBdT = opac.DThermalDistributionOfTNuDT(T, nu);
          kappaPlanckNum += weight * alpha / rho * B * nu * dlnu;
          kappaPlanckDenom += weight * B * nu * dlnu;

          if (alpha > singularity_opac::robust::SMALL()) {
            kappaRosselandNum += weight *
                                 singularity_opac::robust::ratio(rho, alpha) *
                                 dBdT * nu * dlnu;
            kappaRosselandDenom += weight * dBdT * nu * dlnu;
          }
        }

        Real kappaPlanck =
            singularity_opac::robust::ratio(kappaPlanckNum, kappaPlanckDenom);
        Real kappaRosseland = kappaPlanck > singularity_opac::robust::SMALL()
                                  ? singularity_opac::robust::ratio(
                                        kappaRosselandDenom, kappaRosselandNum)
                                  : 0.;
        lkappaPlanck_(iRho, iT) = toLog_(kappaPlanck);
        lkappaRosseland_(iRho, iT) = toLog_(kappaRosseland);
        if (std::isnan(lkappaPlanck_(iRho, iT)) ||
            std::isnan(lkappaRosseland_(iRho, iT))) {
          OPAC_ERROR("photons::MeanOpacity: NAN in opacity evaluations");
        }
      }
    }
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
    lkappaRosseland_.resize(NRho, NT);
    lkappaPlanck_.resize(NRho, NT);

    // set rho-T rankges (Zhu tables are uniform in log-log rho-T space)
    const Real lTMin = toLog_(1.0);
    const Real lTMax = toLog_(7943282.347242886);
    const Real lRhoMin = toLog_(1.0e-14);
    const Real lRhoMax = toLog_(0.7943282347241912);
    lkappaRosseland_.setRange(1, lTMin, lTMax, NT);
    lkappaRosseland_.setRange(2, lRhoMin, lRhoMax, NRho);
    lkappaPlanck_.setRange(1, lTMin, lTMax, NT);
    lkappaPlanck_.setRange(2, lRhoMin, lRhoMax, NRho);

    // fill tables
    for (int iRho = 0; iRho < NRho; ++iRho) {
      const Real lRho_i = lkappaRosseland_.range(2).x(iRho);
      for (int iT = 0; iT < NT; ++iT) {

        // get new file-line
        std::getline(ff, fline);
        cfline = const_cast<char*>(fline.c_str());
        fl_tok = std::strtok(cfline, " ");

        // check for consistent density [g/cm^3] on table row
        const Real Rho = std::stod(fl_tok);
        if (std::abs(Rho - fromLog_(lRho_i)) > 1e-6 * std::abs(Rho)) {
          OPAC_ERROR("photons::MeanOpacity: invalid rho");
        }

        // check for consistent temperature [K] on table row
        const Real lT_i = lkappaRosseland_.range(1).x(iT);
        fl_tok = std::strtok(nullptr, " ");
        const Real T = std::stod(fl_tok);
        if (std::abs(T - fromLog_(lT_i)) > 1e-6 * std::abs(T)) {
          OPAC_ERROR("photons::MeanOpacity: invalid T");
        }

        // populate Rosseland opacity [cm^2/g]
        fl_tok = std::strtok(nullptr, " ");
        lkappaRosseland_(iRho, iT) = toLog_(std::stod(fl_tok));

        // populate Planck opacity [cm^2/g]
        fl_tok = std::strtok(nullptr, " ");
        lkappaPlanck_(iRho, iT) = toLog_(std::stod(fl_tok));

        if (std::isnan(lkappaPlanck_(iRho, iT)) ||
            std::isnan(lkappaRosseland_(iRho, iT))) {
          OPAC_ERROR("photons::MeanOpacity: NAN in parsed ASCII opacity");
        }
      }
    }
  }

  PORTABLE_INLINE_FUNCTION Real toLog_(const Real x) const {
    return std::log10(std::abs(x) + EPS);
  }
  PORTABLE_INLINE_FUNCTION Real fromLog_(const Real lx) const {
    return std::pow(10., lx);
  }
  Spiner::DataBox<Real> lkappaPlanck_;
  Spiner::DataBox<Real> lkappaRosseland_;
  const char *filename_;
};

#undef EPS

} // namespace impl

using MeanOpacityBase = impl::MeanOpacity;
using MeanOpacity =
    impl::MeanVariant<MeanOpacityBase, MeanNonCGSUnits<MeanOpacityBase>>;

} // namespace photons
} // namespace singularity

#endif // SINGULARITY_OPAC_PHOTONS_MEAN_OPACITY_PHOTONS__
