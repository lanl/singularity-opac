// ======================================================================
// © 2026. Triad National Security, LLC. All rights reserved.  This
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
#ifndef SINGULARITY_OPAC_BASE_SP5_
#define SINGULARITY_OPAC_BASE_SP5_
// This file was made in part with generative AI.

namespace SP5 {

namespace Opac {
constexpr char defaultFileName[] = "opac.sp5";
constexpr char AbsorptionCoefficient[] = "absorption coefficient";
constexpr char AngleAveragedAbsorptionCoefficient[] =
    "angle-averaged absorption coefficient";
constexpr char EmissivityPerNu[] = "emissivity per nu";
constexpr char TotalEmissivity[] = "total emissivity";
constexpr char NumberEmissivity[] = "number emissivity";
} // namespace Opac

namespace MeanOpac {
constexpr char PlanckMeanOpacity[] = "Planck mean opacity";
constexpr char RosselandMeanOpacity[] = "Rosseland mean opacity";
} // namespace MeanOpac

namespace MeanSOpac {
constexpr char PlanckMeanSOpacity[] = "Planck mean scattering opacity";
constexpr char RosselandMeanSOpacity[] = "Rosseland mean scattering opacity";
} // namespace MeanSOpac

namespace Multigroup {
constexpr char GroupBounds[] = "group bounds";
} // namespace Multigroup

namespace MultigroupOpac {
constexpr char PlanckGroupOpacity[] = "Planck group opacity";
constexpr char RosselandGroupOpacity[] = "Rosseland group opacity";
} // namespace MultigroupOpac

namespace MultigroupSOpac {
constexpr char PlanckGroupSOpacity[] = "Planck group scattering opacity";
constexpr char RosselandGroupSOpacity[] = "Rosseland group scattering opacity";
} // namespace MultigroupSOpac

namespace Offsets {
constexpr char opac_messageName[] = "opac_interpretation";
constexpr char opac_message[] =
    "Opacity quantities are functions of log_10(X)\n"
    "for X = density rho or temperature T; group boundaries nu are physical\n"
    "frequency values in Hz and are indexed directly.\n";
constexpr char opac_rho[] = "opac_rhoOffset";
constexpr char opac_T[] = "opac_TOffset";
} // namespace Offsets

namespace Material {
constexpr char opac_comments[] = "opac_comments";
constexpr char opacid[] = "opacid";
constexpr char opac_name[] = "opac_name";
} // namespace Material

// IPCRESS-only quantities
namespace IPCRESS {
constexpr char RosselandTotalMultigroupOpacity[] =
    "rosseland total multigroup opacity";
constexpr char PlanckTotalMultigroupOpacity[] =
    "planck total multigroup opacity";
constexpr char RosselandTotalGrayOpacity[] = "rosseland total gray opacity";
constexpr char PlanckTotalGrayOpacity[] = "planck total gray opacity";
} // namespace IPCRESS

} // namespace SP5

#endif // SINGULARITY_OPAC_BASE_SP5_
