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
(here, $\sigma$ is the frequency- and angle-dependent cross section in units of $cm^2$):
| Function              | Expression | Description            | Units   |
| --------------------- | ---------- | ---------------------  | ------- |
| AbsorptionCoefficient | $n \sigma$ | Absorption coefficient | $cm^-1$ |


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
| SINGULARITY_USE_HDF5              | ON      | Enables HDF5. Required for SpinerEOS, StellarCollapseEOS, and sesame2spiner          |
| SINGULARITY_USE_FORTRAN           | ON      | Enable Fortran API for equation of state                                             |
| SINGULARITY_USE_KOKKOS            | OFF     | Uses Kokkos as the portability backend. Currently only Kokkos is supported for GPUs. |
| SINGULARITY_USE_EOSPAC            | OFF     | Link against EOSPAC. Needed for sesame2spiner and some tests.                        |
| SINGULARITY_USE_CUDA              | OFF     | Target nvidia GPUs via cuda. Currently requires Kokkos.                              |
| SINGULARITY_USE_KOKKOSKERNELS     | OFF     | Use Kokkos Kernels for linear algebra. Needed for mixed cell closure models on GPU   |
| SINGULARITY_BUILD_CLOSURE         | ON      | Builds mixed cell closure machinery for multi-material problems                      |
| SINGULARITY_BUILD_TESTS           | OFF     | Build test infrastructure.                                                           |
| SINGULARITY_TEST_SESAME           | OFF     | Test the Sesame table readers                                                        |
| SINGULARITY_TEST_STELLAR_COLLAPSE | OFF     | Test the Stellar Collapse table readers                                              |
| SINGULARITY_BUILD_SESAME2SPINER   | OFF     | Builds the conversion tool sesame2spiner which makes files readable by SpinerEOS     |
| SINGULARITY_BUILD_STELLARCOLLAPSE2SPINER | OFF     | Builds the conversion tool stellarcollapse2spiner which optionally makes stellar collapse files faster to read |
| SINGULARITY_INVERT_AT_SETUP       | OFF     | For tests, pre-invert eospac tables.                                                 |
| SINGULARITY_BETTER_DEBUG_FLAGS    | ON      | Enables nicer GPU debug flags. May interfere with in-tree builds as a submodule      |
| SINGULARITY_HIDE_MORE_WARNINGS    | OFF     | Makes warnings less verbose. May interfere with in-tree builds as a submodule        |
| SINGULARITY_SUBMODULE_MODE        | OFF     | Sets cmake flags for in-tree builds                                                  |

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
