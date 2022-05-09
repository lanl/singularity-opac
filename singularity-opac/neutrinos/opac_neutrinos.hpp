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

#ifndef SINGULARITY_OPAC_NEUTRINOS_OPAC_NEUTRINOS_
#define SINGULARITY_OPAC_NEUTRINOS_OPAC_NEUTRINOS_

#include <variant/include/mpark/variant.hpp>

#include <singularity-opac/neutrinos/brt_neutrinos.hpp>
#include <singularity-opac/neutrinos/gray_opacity_neutrinos.hpp>
#include <singularity-opac/neutrinos/neutrino_variant.hpp>
#include <singularity-opac/neutrinos/non_cgs_neutrinos.hpp>
#include <singularity-opac/neutrinos/spiner_opac_neutrinos.hpp>
#include <singularity-opac/neutrinos/thermal_distributions_neutrinos.hpp>
#include <singularity-opac/neutrinos/tophat_emissivity_neutrinos.hpp>

#include <singularity-opac/neutrinos/mean_opacity_neutrinos.hpp>

namespace singularity {
namespace neutrinos {

// TODO(JMM): Include chemical potential
using ScaleFree =
    GrayOpacity<FermiDiracDistributionNoMu<3, PhysicalConstantsUnity>,
                PhysicalConstantsUnity>;
using BRTOpac = BRTOpacity<FermiDiracDistributionNoMu<3>>;
using Gray = GrayOpacity<FermiDiracDistributionNoMu<3>>;
using Tophat = TophatEmissivity<FermiDiracDistributionNoMu<3>>;
using SpinerOpac = SpinerOpacity<FermiDiracDistributionNoMu<3>>;

using Opacity = impl::Variant<ScaleFree, BRTOpac, Gray, Tophat, SpinerOpac,
                              NonCGSUnits<BRTOpac>, NonCGSUnits<Gray>,
                              NonCGSUnits<Tophat>, NonCGSUnits<SpinerOpac>>;

} // namespace neutrinos
} // namespace singularity

#endif // SINGULARITY_OPAC_NEUTRINOS_OPAC_NEUTRINOS_
