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

namespace MultigroupOpac {
constexpr char PlanckGroupOpacity[] = "Planck group opacity";
constexpr char RosselandGroupOpacity[] = "Rosseland group opacity";
constexpr char GroupBounds[] = "group bounds";
} // namespace MultigroupOpac

namespace MultigroupSOpac {
constexpr char PlanckGroupSOpacity[] = "Planck group scattering opacity";
constexpr char RosselandGroupSOpacity[] = "Rosseland group scattering opacity";
constexpr char GroupBounds[] = "group bounds";
} // namespace MultigroupSOpac

// constants and fields below used in ipcress2spiner utility
constexpr char defaultSesFileName[] = "converted_ipcress.sp5";
constexpr char logType[] = "log_type";

namespace Depends {
constexpr char logRhoLogSie[] = "dependsLogRhoLogSie";
constexpr char logRhoLogT[] = "dependsLogRhoLogT";
} // namespace Depends

namespace SubTable {
constexpr char electronOnly[] = "electronOnly";
constexpr char ionCold[] = "ionCold";
} // namespace SubTable

namespace Offsets {
constexpr char messageName[] = "interpretation";
constexpr char message[] =
    "All quantities are functions of log_10(X)\n"
    "for X = density rho, temperature T, or group boundaries hnu\n"
    "where conversion is X = 10^{Xlog}\n";
constexpr char rho[] = "rhoOffset";
constexpr char T[] = "TOffset";
constexpr char group_bounds[] = "groupBoundsOffset";
} // namespace Offsets

namespace Material {
constexpr char comments[] = "comments";
constexpr char matid[] = "matid";
constexpr char name[] = "name";
} // namespace Material

namespace Fields {
constexpr char ramg[] = "rosseland absorption multigroup opacity";
constexpr char rsmg[] =  "rosseland scattering multigroup opacity";
constexpr char rtmg[] =  "rosseland total multigroup opacity";
constexpr char pmg[] = "planck total multigroup opacity";
constexpr char ragray[] = "rosseland absorption gray opacity";
constexpr char rgray[] = "rosseland total gray opacity";
constexpr char pgray[] =  "planck total gray opacity";
constexpr char T[] = "temperature";
constexpr char rho[] = "density";
constexpr char group_bounds[] = "group boundaries";
} // namespace Fields

} // namespace SP5

#endif // SINGULARITY_OPAC_BASE_SP5_
