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

#ifndef SINGULARITY_OPAC_PHOTONS_OPAC_PHOTONS_
#define SINGULARITY_OPAC_PHOTONS_OPAC_PHOTONS_

#include <singularity-opac/photons/epbremsstrahlung_opacity_photons.hpp>
#include <singularity-opac/photons/gray_opacity_photons.hpp>
#include <singularity-opac/photons/non_cgs_photons.hpp>
#include <singularity-opac/photons/photon_variant.hpp>
#include <singularity-opac/photons/powerlaw_opacity_photons.hpp>
#include <singularity-opac/photons/thermal_distributions_photons.hpp>

#include <singularity-opac/photons/mean_opacity_photons.hpp>

namespace singularity {
namespace photons {

using ScaleFree = GrayOpacity<PhysicalConstantsUnity>;
using Gray = GrayOpacity<PhysicalConstantsCGS>;
using PowerLawScaleFree = PowerLawOpacity<PhysicalConstantsUnity>;
using PowerLaw = PowerLawOpacity<PhysicalConstantsCGS>;
using EPBremss = EPBremsstrahlungOpacity<PhysicalConstantsCGS>;

using Opacity = impl::Variant<ScaleFree, Gray, PowerLawScaleFree, PowerLaw,
                              EPBremss, NonCGSUnits<Gray>,
                              NonCGSUnits<PowerLaw>, NonCGSUnits<EPBremss>>;

} // namespace photons
} // namespace singularity

#endif // SINGULARITY_OPAC_PHOTONS_OPAC_PHOTONS_
