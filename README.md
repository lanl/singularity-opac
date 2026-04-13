<!-- This file was made in part with generative AI. -->
singularity-opac
==

[![Tests](https://github.com/lanl/singularity-opac/actions/workflows/tests.yml/badge.svg)](https://github.com/lanl/singularity-opac/actions/workflows/tests.yml)

Performance Portable Opacity and Emissivity library for simulation codes

## API

singularity-opac provides a uniform API for all opacity models, in two forms: frequency-dependent, and frequency-averaged (Plank or Rosseland means), and separately for absorption and scattering opacities.

For frequency-dependent absorption opacities, the following functions are provided
(here, $\sigma$ is the frequency- and angle-dependent cross section in units of ${\rm cm}^2$):
| Function              | Expression | Description            | Units   |
| --------------------- | ---------- | ---------------------  | ------- |
| AbsorptionCoefficient | $n \sigma$ | Absorption coefficient | ${\rm cm}^{-1}$ |
| AngleAveragedAbsorptionCoefficient | $\frac{1}{4 \pi}\int n \sigma d\Omega$ | Absorption coefficient averaged over solid angle | ${\rm cm}^{-1}$ |
| EmissivityPerNuOmega | $j_{\nu} = \frac{dE}{d^3x dt d\Omega d\nu}$ | Frequency- and angle-dependent emissivity | ${\rm erg}{\rm cm}^{-3}{\rm s}^{-1}{\rm Hz}^{-1}{\rm Sr}^{-1}$ |
| EmissivityPerNu | $\int j_{\nu} d\Omega$  | Frequency-dependent emissivity | ${\rm erg}{\rm cm}^{-3}{\rm s}^{-1}{\rm Hz}^{-1}$ |
| Emissivity | $\int j_{\nu} d\nu d\Omega$  | Total emissivity | ${\rm erg}{\rm cm}^{-3}{\rm s}^{-1}$ |
| NumberEmissivity | $\int \frac{1}{h \nu} j_{\nu} d\Omega d\nu$ | Total number emissivity | ${\rm cm}^{-3}{\rm s}^{-1}$ |
| ThermalDistributionOfTNu | $B_{\nu} = \frac{dE}{dA dt d\Omega d\nu}$ | Specific intensity of thermal distribution | ${\rm erg}{\rm cm}^{-2}{\rm s}^{-1}{\rm Sr}^{-1}{\rm Hz}^{-1}$ |
| DThermalDistributionOfTNuDT | $dB_{\nu}/dT$ | Temperature derivative of specific intensity of thermal distribution | ${\rm erg}{\rm cm}^{-2}{\rm s}^{-1}{\rm Sr}^{-1}{\rm Hz}^{-1}{\rm K}^{-1}$ |
| ThermalDistributionOfT | $4 \pi \int B_{\nu} d\nu = \int B_{\nu} d\Omega d\nu$ | Frequency- and angle-integrated intensity of thermal distribution | ${\rm erg}{\rm cm}^{-2}{\rm s}^{-1}$ |
| ThermalNumberDistributionOfT | $B = \int \frac{1}{h \nu} B_{\nu} d\Omega d\nu$ | Frequency- and angle-integrated intensity of thermal distribution | ${\rm erg}{\rm cm}^{-2}{\rm s}^{-1}$ |
| EnergyDensityFromTemperature | $E_{\rm R}$ | Radiation energy density | ${\rm erg}{\rm cm}^{-3}$ |
| TemperatureFromEnergyDensity | $T_{\rm R}$ | Radiation temperature | ${\rm K}$ |
| NumberDensityFromTemperature | $n_{\rm R}$ | Radiation number density | ${\rm cm}^{-3}$ |

with the following function signatures:

    AbsorptionCoefficient(density, temperature, frequency)
    AngleAveragedAbsorptionCoefficient(density, temperature, frequency)
    EmissivityPerNuOmega(density, temperature, frequency)
    EmissivityPerNu(density, temperature, frequency)
    Emissivity(density, temperature)
    NumberEmissivity(density, temperature)
    ThermalDistributionOfTNu(temperature, frequency)
    DThermalDistribtuionOfTNuDT(temperature, frequency)
    ThermalDistributionOfT(temperature)
    ThermalNumberDistributionOfT(temperature)
    EnergyDensityFromTemperature(temperature)
    TemperatureFromEnergyDensity(radiation energy density)
    NumberDensityFromTemperature(temperature)

For mean absorption opacities, the following functions are provided.
Mean opacities compute frequency-averaged (Planck or Rosseland) absorption
coefficients over one or more energy groups. When `ngroups=1` with group bounds
`[0, ∞)`, mean opacities reduce to traditional gray opacities. For `ngroups>1`,
they provide multigroup radiation transport capabilities.

| Function              | Expression | Description            | Units   |
| --------------------- | ---------- | ---------------------  | ------- |
| PlankMeanAbsorptionCoefficient | $n \sigma$ | Planck mean absorption coefficient (gray, ngroups=1) | ${\rm cm}^{-1}$ |
| RosselandMeanAbsorptionCoefficient | $n \sigma$ | Rosseland mean absorption coefficient (gray, ngroups=1) | ${\rm cm}^{-1}$ |
| PlanckGroupAbsorptionCoefficient | $n \sigma_g$ | Planck-weighted absorption coefficient in a frequency group | ${\rm cm}^{-1}$ |
| RosselandGroupAbsorptionCoefficient | $n \sigma_g$ | Rosseland-weighted absorption coefficient in a frequency group | ${\rm cm}^{-1}$ |
| AbsorptionCoefficient | $n \sigma$ or $n \sigma_g$ | Absorption coefficient (gray or group) | ${\rm cm}^{-1}$ |
| Emissivity | $\int j_{\nu} d\nu d\Omega$  | Total emissivity | ${\rm erg}{\rm cm}^{-3}{\rm s}^{-1}$ |
| GroupOfNu | $g(\nu)$ | Group index containing a given frequency | dimensionless |
| PlanckGroupAbsorptionCoefficientFromNu | $n \sigma_{g(\nu)}$ | Planck-weighted group coefficient selected by frequency | ${\rm cm}^{-1}$ |
| RosselandGroupAbsorptionCoefficientFromNu | $n \sigma_{g(\nu)}$ | Rosseland-weighted group coefficient selected by frequency | ${\rm cm}^{-1}$ |
| AbsorptionCoefficientFromNu | $n \sigma_{g(\nu)}$ | Group coefficient selected by frequency | ${\rm cm}^{-1}$ |

with the following function signatures:

    ngroups()
    PlanckMeanAbsorptionCoefficient(density, temperature)
    RosselandMeanAbsorptionCoefficient(density, temperature)
    PlanckGroupAbsorptionCoefficient(density, temperature, group index)
    RosselandGroupAbsorptionCoefficient(density, temperature, group index)
    AbsorptionCoefficient(density, temperature, gmode [Planck, Rosseland])
    AbsorptionCoefficient(density, temperature, group index, gmode [Planck, Rosseland])
    Emissivity(density, temperature)
    GroupOfNu(frequency)
    PlanckGroupAbsorptionCoefficientFromNu(density, temperature, frequency)
    RosselandGroupAbsorptionCoefficientFromNu(density, temperature, frequency)
    AbsorptionCoefficientFromNu(density, temperature, frequency, gmode [Planck, Rosseland])

Mean absorption opacities can either be constructed from a frequency-dependent
Mean absorption opacities can either be constructed from a frequency-dependent
opacity model plus user-provided group bounds, or loaded directly from
precomputed Spiner tables of group Planck and Rosseland opacities. The direct
table-backed path does not require singularity-opac to recompute opacity
integrals, but it must include explicit group bounds. Mean opacities
always carry group bounds, and the convention is half-open groups
`[nu_g, nu_{g+1})`, with the final upper bound included in the last group.
If you want singularity-opac to interpret a group as extending to infinity,
then the final entry of the `group bounds` array must literally be IEEE
positive infinity, not just a very large finite cutoff. In C++, that means
using `std::numeric_limits<Real>::infinity()`, while in Python that would
typically be `float("inf")` or `numpy.inf`. In an HDF/Spiner table, that same
IEEE `+infinity` value must appear in the final element of the `"group bounds"`
dataset. Likewise, a lower tail group `[0, nu_1)` is represented by setting
the first bound to `0.`. A very large finite number is still interpreted as a
finite bound.

For frequency-dependent scattering opacities, the following functions are provided
| Function              | Expression | Description            | Units   |
| --------------------- | ---------- | ---------------------  | ------- |
| TotalCrossSection | $\sigma$ | Scattering cross section | ${\rm cm}^{2}$ |
| DifferentialCrossSection | $d\sigma / d \Omega $ | Differential scattering cross section | ${\rm cm}^{2}{\rm Sr}^{-1}$ |
| TotalScatteringCoefficient | $n \sigma $ | Scattering coefficient | ${\rm cm}^{-1}$ |

with the following function signatures:

    TotalCrossSection(density, temperature, frequency)
    DifferentialCrossSection(density, temperature, frequency, cos(theta))
    TotalScatteringCoefficient(density, temperature, frequency)

For mean scattering opacities, the following functions are provided:
| Function              | Expression | Description            | Units   |
| --------------------- | ---------- | ---------------------  | ------- |
| PlanckMeanScatteringCoefficient | $n \sigma$ | Planck mean scattering coefficient | ${\rm cm}^{-1}$ |
| RosselandMeanScatteringCoefficient | $n \sigma$ | Rosseland mean scattering coefficient | ${\rm cm}^{-1}$ |

with the following function signatures:

    PlanckMeanScatteringCoefficient(density, temperature)
    RosselandMeanScatteringCoefficient(density, temperature)

Note that `ThermalDistributionOfTNu` is the per-steradian Planck function `B_\nu`, so
`\int B_\nu d\nu = (c / 4 \pi) a T^4`, while
`ThermalDistributionOfT = 4 \pi \int B_\nu d\nu = c a T^4`. Therefore the
thermal radiation energy density is
`u = (4 \pi / c) \int B_\nu d\nu = (1 / c) ThermalDistributionOfT = a T^4`, and the
thermal radiation number density is
`n = (1 / c) ThermalNumberDistributionOfT`.

Opacity variant constructors are specific to the opacity model being requested; consult the source code for
individual opacities.

Internally singularity-opac always uses CGS units, as in the above table. However, arbitrary units are supported through the units modifier, which accepts
function argument inputs in the arbitrary unit system, and returns the result from the function in those same arbitrary units. For example, a gray absorption opacity in non-cgs units specified by `time_unit`, `mass_unit`, `length_unit`, and `temp_unit` conversion factors from code to CGS units (e.g. `mass_cgs = mass_unit * mass_code`)
is created as

    photons::Opacity noncgs_opacity = photons::NonCGSUnits<photons::Gray>(
          photons::Gray(kappa), time_unit, mass_unit, length_unit, temp_unit);

Note that neutrino opacity functions also include electron fraction and RadiationType species arguments.

Frequency-dependent emissition and absorption functions do not currently support angle dependence.

A struct of runtime physical constants is provided for optional consistency with internal operations by the
`GetRuntimePhysicalConstants()` method.

## To Build

At its most basic:
```bash
git clone --recursive git@github.com:lanl/singularity-opac.git
cd singularity-opac
mkdir bin
cd bin
cmake ..
make
```

### To Make tests

Then, in the singularity-opac root directory:
```bash
mkdir -p bin
cd bin
cmake -DSINGULARITY_BUILD_TESTS=ON ..
make -j
make test
```

### Build Options

A number of options are avaialable for compiling:

| Option                            | Default | Comment                                                                              |
| --------------------------------- | ------- | ------------------------------------------------------------------------------------ |
| SINGULARITY_BUILD_TESTS           | OFF     | Build test infrastructure.                                                           |
| SINGULARITY_USE_HDF5              | ON      | Enables HDF5. Required for Spiner opacities.                                         |
| SINGULARITY_KOKKOS_IN_TREE        | OFF     | Force cmake to use Kokkos source included in tree.                                   |

### Loading ASCII Data

Currently, the MeanOpacity class, defined in singularity-opac/photons/mean_opacity_photons.hpp, supports
loading grey Rosseland and Planck opacity data in an ASCII format.  An example of this format is
provided by singularity-opac/photons/example_ascii/kap_plaw.txt.  The 1st row of the header has the
number of density points, NRho, then the number of temperature points, NT.  The 2nd (3rd) row of the header
has min and max density (temperature) bounds.  These bounds are inclusive, so the opacity data in the file
should have evaluations at these min and max values.  The rest of the ASCII file is a two-column table, where
the 1st (2nd) column is Rosseland (Planck) opacity.  The number of rows in each column is NRhoxNT, where
density is the slow index and temperature is the fast index (thus the row index = temperature index
+ NT x (density index), indexing from 0).  Each opacity is assumed to be evaluated on log-spaced
density and temperature grids, where these grids are defined by NRho, NT, and the (again inclusive) min and
max bounds the header.

## Copyright

© 2021-2026. Triad National Security, LLC. All rights reserved.  This
program was produced under U.S. Government contract 89233218CNA000001
for Los Alamos National Laboratory (LANL), which is operated by Triad
National Security, LLC for the U.S.  Department of Energy/National
Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of
Energy/National Nuclear Security Administration. The Government is
granted for itself and others acting on its behalf a nonexclusive,
paid-up, irrevocable worldwide license in this material to reproduce,
prepare derivative works, distribute copies to the public, perform
publicly and display publicly, and to permit others to do so.
