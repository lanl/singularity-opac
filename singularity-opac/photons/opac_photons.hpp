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

#include <singularity-opac/base/opac_variant.hpp>

#include <singularity-opac/photons/gray_opacity_photons.hpp>
#include <singularity-opac/photons/thermal_distributions_photons.hpp>

namespace singularity {
namespace photons {

using Gray = GrayOpacity;
using Opacity = opac_impl::Variant<Gray>;

} // namespace photons
} // namespace singularity

#endif // SINGULARITY_OPAC_PHOTONS_OPAC_PHOTONS_
