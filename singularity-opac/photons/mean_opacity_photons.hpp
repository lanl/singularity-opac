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

#ifndef SINGULARITY_OPAC_PHOTONS_MEAN_OPACITY_PHOTONS_
#define SINGULARITY_OPAC_PHOTONS_MEAN_OPACITY_PHOTONS_

#include <cmath>

#include <ports-of-call/portability.hpp>
#include <singularity-opac/base/radiation_types.hpp>
#include <singularity-opac/base/sp5.hpp>
#include <singularity-opac/constants/constants.hpp>
#include <spiner/databox.hpp>

#include <singularity-opac/photons/mean_photon_variant.hpp>

namespace singularity {
namespace photons {
namespace impl {

#define EPS (10.0 * std::numeric_limits<Real>::min())

// TODO(BRR) Note: It is assumed that lambda is constant for all densities,
// temperatures, and Ye

template <typename pc = PhysicalConstantsCGS>
class MeanOpacity {
 public:
  MeanOpacity() = default;
  template <typename Opacity>
  MeanOpacity(const Opacity &opac, const Real lRhoMin, const Real lRhoMax,
              const int NRho, const Real lTMin, const Real lTMax, const int NT,
              Real *lambda = nullptr) {
    lkappaPlanck_.resize(NRho, NT);
    // index 0 is the species and is not interpolatable
    lkappaPlanck_.setRange(1, lTMin, lTMax, NT);
    lkappaPlanck_.setRange(2, lRhoMin, lRhoMax, NRho);
    lkappaRosseland_.copyMetadata(lkappaPlanck_);

    // Fill tables
    for (int iRho = 0; iRho < NRho; ++iRho) {
      Real lRho = lkappaPlanck_.range(3).x(iRho);
      Real rho = fromLog_(lRho);
      for (int iT = 0; iT < NT; ++iT) {
        Real lT = lkappaPlanck_.range(2).x(iT);
        Real T = fromLog_(lT);
            Real kappaPlanckNum = 0.;
            Real kappaPlanckDenom = 0.;
            Real kappaRosselandNum = 0.;
            Real kappaRosselandDenom = 0.;
            // Integrate over frequency
            const int nnu = 100;
            const Real lnuMin =
                toLog_(1.e-3 * pc::kb * fromLog_(lTMin) / pc::h);
            const Real lnuMax = toLog_(1.e3 * pc::kb * fromLog_(lTMax) / pc::h);
            const Real dlnu = (lnuMax - lnuMin) / (nnu - 1);
            for (int inu = 0; inu < nnu; ++inu) {
              const Real lnu = lnuMin + inu * dlnu;
              const Real nu = fromLog_(lnu);
              kappaPlanckNum +=
                  opac.AbsorptionCoefficient(rho, T, nu, lambda) /
                  rho * opac.ThermalDistributionOfTNu(T, nu) * nu * dlnu;
              kappaPlanckDenom +=
                  opac.ThermalDistributionOfTNu(T, nu) * nu * dlnu;

              kappaRosselandNum +=
                  rho /
                  opac.AbsorptionCoefficient(rho, T, nu, lambda) *
                  opac.DThermalDistributionOfTNuDT(T, nu) * nu * dlnu;
              kappaRosselandDenom +=
                  opac.DThermalDistributionOfTNuDT(T, nu) * nu * dlnu;
            }

            // Trapezoidal rule
            const Real nu0 = fromLog_(lnuMin);
            const Real nu1 = fromLog_(lnuMax);
            kappaPlanckNum -=
                0.5 * 1. / rho *
                (opac.AbsorptionCoefficient(rho, T, nu0, lambda) *
                     opac.ThermalDistributionOfTNu(T, nu0) * nu0 +
                 opac.AbsorptionCoefficient(rho, T, nu1, lambda) *
                     opac.ThermalDistributionOfTNu(T, nu1) * nu1) *
                dlnu;
            kappaPlanckDenom -=
                0.5 *
                (opac.ThermalDistributionOfTNu(T, nu0) * nu0 +
                 opac.ThermalDistributionOfTNu(T, nu1) * nu1) *
                dlnu;
            kappaRosselandNum -=
                0.5 * rho *
                (1. /
                     opac.AbsorptionCoefficient(rho, T, nu0, lambda) *
                     opac.DThermalDistributionOfTNuDT(T, nu0) * nu0 +
                 1. /
                     opac.AbsorptionCoefficient(rho, T, nu1, lambda) *
                     opac.DThermalDistributionOfTNuDT(T, nu1) * nu1) *
                dlnu;
            kappaRosselandDenom -=
                0.5 *
                (opac.DThermalDistributionOfTNuDT(T, nu0) * nu0 +
                 opac.DThermalDistributionOfTNuDT(T, nu1) * nu1) *
                dlnu;

            Real lkappaPlanck = toLog_(kappaPlanckNum / kappaPlanckDenom);
            Real lkappaRosseland =
                toLog_(1. / (kappaRosselandNum / kappaRosselandDenom));
            lkappaPlanck_(iRho, iT) = lkappaPlanck;
            lkappaRosseland_(iRho, iT) = lkappaRosseland;
            if (std::isnan(lkappaPlanck_(iRho, iT)) ||
                std::isnan(lkappaRosseland_(iRho, iT))) {
              OPAC_ERROR("photons::MeanOpacity: NAN in opacity evaluations");
            }
          }
        }
      }
    }
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
  Real RosselandMeanAbsorptionCoefficient(const Real rho, const Real temp) const {
    Real lRho = toLog_(rho);
    Real lT = toLog_(temp);
    return rho * fromLog_(lkappaRosseland_.interpToReal(lRho, lT));
  }

 private:
  PORTABLE_INLINE_FUNCTION Real toLog_(const Real x) const {
    return std::log10(std::abs(x) + EPS);
  }
  PORTABLE_INLINE_FUNCTION Real fromLog_(const Real lx) const {
    return std::pow(10., lx);
  }
  Spiner::DataBox lkappaPlanck_;
  Spiner::DataBox lkappaRosseland_;
  const char *filename_;
};

#undef EPS

} // namespace impl

using MeanOpacityScaleFree = impl::MeanOpacity<PhysicalConstantsUnity>;
using MeanOpacityCGS = impl::MeanOpacity<PhysicalConstantsCGS>;
using MeanOpacity = impl::MeanVariant<MeanOpacityScaleFree, MeanOpacityCGS>;

} // namespace photons
} // namespace singularity

#endif // SINGULARITY_OPAC_PHOTONS_MEAN_OPACITY_PHOTONS__
