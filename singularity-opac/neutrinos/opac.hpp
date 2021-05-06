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

#ifndef OPACITIES_OPAC_
#define OPACITIES_OPAC_

#include <singularity-opac/base/opac_variant.hpp>

#include "gray_opacity.hpp"
#include "thermal_distributions.hpp"
#include "tophat_emissivity.hpp"

namespace singularity {

using GrayPhotonOpacity = GrayOpacity<PlanckDistribution<1>>;
// TODO(JMM): Change this to Fermi-Dirac
using GrayNeutrinoOpacity = GrayOpacity<PlanckDistribution<3>>;
using TophatOpacity = TophatEmissivity<PlanckDistribution<3>>;

using Opacity =
    opac_impl::Variant<GrayPhotonOpacity, GrayNeutrinoOpacity, TophatOpacity>;

} // namespace singularity

#endif // OPACITIES_OPAC_
