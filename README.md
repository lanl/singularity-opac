singularity-opac
==

[![Tests](https://github.com/lanl/singularity-opac/actions/workflows/tests.yml/badge.svg)](https://github.com/lanl/singularity-opac/actions/workflows/tests.yml)

Performance Portable Opacity and Emissivity for simulation codes

## Components

- `eos` contains the analytic and tabular equation of state infrastructure
- `eospac-wrapper` is a wrapper for eospac
- `closure` contains mixed cell closure models. Currently pressure-temperature equilibrium is supported. Requires `eos`.
- `sesame2spiner` converts sesame tables to spiner tables
- `stellarcollaspe2spienr` converts tables found on stellarcollapse.org to spiner tables

## API

singularity-opac provides a uniform API for all opacity models. The following functions are provided
(here, `\sigma` is the frequency- and angle-dependent cross section in units of `cm^2`):
| Function              | Expression | Description            | Units   |
| --------------------- | ---------- | ---------------------  | ------- |
| AbsorptionCoefficient | `n \sigma` | Absorption coefficient | `cm^{-1}` |
| AngleAveragedAbsorptionCoefficient | `\frac{1}{4 \pi}\int d\Omega n \sigma` |Absorption coefficient averaged over solid angle | `cm^{-1}` |
| EmissivityPerNuOmega | `j_{\nu} = \frac{dE}{d^3x dt d\Omega d\nu}` | Frequency- and angle-dependent emissivity | $erg cm^{-3} s^{-1} Sr^{-1} Hz^{-1}$ |
| EmissivityPerNu | `\int j_{\nu} d\nu`  | Frequency-dependent emissivity | `erg cm^{-3} s^{-1} Sr^{-1}` |
| Emissivity | `\int j_{\nu} d\nu d\Omega`  | Total emissivity | `erg cm^{-3} s^{-1}` |

Note that the thermal radiation energy density `u = 1/c ThermalDistributionOfT` and the thermal radiation number density `n = 1/c ThermalNumberDistribution`.


## To Build

At its most basic:
```bash
git clone --recursive git@github.com:lanl/singularity-eos.git
cd singularity-eos
mkdir bin
cd bin
cmake ..
make
```

### To Make tests

Then, in the singularity-eos root directory:
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

## Units

singularity-eos, prior to unit modifiers, always uses CGS units internally. For
concreteness, the internal units of each quantity are:

| Quantity                  | Alternate names      | Units         |
| ------------------------- | -------------------- | ------------- |
| Mass density              | rho                  | g cm^-3       |
| Temperature               | temperature          | K             |
| Specific Internal Energy  | sie, Internal Energy | erg g^-1      |
| Pressure                  |                      | erg cm^-3     |
| Specific Heat             |                      | erg g^-1 K^-1 |
| Bulk Modulus              |                      | erg cm^-3     |
| Gruneisen Parameter       |                      |               |

## Copyright

© 2021. Triad National Security, LLC. All rights reserved.  This
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
