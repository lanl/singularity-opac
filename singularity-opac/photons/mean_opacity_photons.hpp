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

#ifdef SPINER_USE_HDF
  MeanOpacity(const std::string &filename) : filename_(filename.c_str()) {
    herr_t status = H5_SUCCESS;
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    status += lkappaPlanck_.loadHDF(file, SP5::MeanOpac::PlanckMeanOpacity);
    status +=
        lkappaRosseland_.loadHDF(file, SP5::MeanOpac::RosselandMeanOpacity);
    status += H5Fclose(file);

    if (status != H5_SUCCESS) {
      OPAC_ERROR("photons::MeanOpacity: HDF5 error\n");
    }
  }

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

using MeanOpacity =
    impl::MeanVariant<impl::MeanOpacity, MeanNonCGSUnits<impl::MeanOpacity>>;

} // namespace photons
} // namespace singularity

#endif // SINGULARITY_OPAC_PHOTONS_MEAN_OPACITY_PHOTONS__
