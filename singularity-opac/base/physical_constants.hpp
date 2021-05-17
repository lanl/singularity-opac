// ======================================================================
// Copyright © 2018, University of Illinois. All rights reserved.
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

#ifndef SINGULARITY_OPAC_BASE_PHYSICAL_CONSTANTS_
#define SINGULARITY_OPAC_BASE_PHYSICAL_CONSTANTS_

#include <cmath>
#include <ports-of-call/portability.hpp>

// Convention taken from ebhlight/nubhlight
namespace singularity {
namespace constants {
namespace cgs {

constexpr Real EE = 4.80320680e-10;       // Electron charge
constexpr Real CL = 2.99792458e10;        // Speed of light
constexpr Real ME = 9.1093826e-28;        // Electron mass
constexpr Real MP = 1.67262171e-24;       // Proton mass
constexpr Real MN = 1.67492728e-24;       // Neutron mass
constexpr Real HPL = 6.6260693e-27;       // Planck constant
constexpr Real HBAR = HPL / (2. * M_PI);  // Reduced Planck constant
constexpr Real KBOL = 1.3806505e-16;      // Boltzmann constant
constexpr Real GNEWT = 6.6742e-8;         // Gravitational constant
constexpr Real SIG = 5.670400e-5;         // Stefan-Boltzmann constant
constexpr Real AR = 4 * SIG / CL;         // Radiation constant
constexpr Real THOMSON = 0.665245873e-24; // Thomson cross section
constexpr Real COULOMB_LOG = 20.;         // Coulomb logarithm
constexpr Real ALPHAFS = 0.007299270073;  // Fine structure constant ~ 1./137.
constexpr Real GFERM = 1.435850814e-49;   // Fermi constant
constexpr Real GA = -1.272323;            // Axial-vector coupling
constexpr Real GA2 = GA * GA;
constexpr Real S2THW = 0.222321; // sin^2(Theta_W), Theta_W = Weinberg angle
constexpr Real S4THW = S2THW * S2THW;
constexpr Real NUSIGMA0 =
    1.7611737037e-44; // Fundamental neutrino cross section

// Unit conversions
constexpr Real EV = 1.60217653e-12;   // Electron-volt
constexpr Real MEV = 1.0e6 * EV;      // Mega-Electron-Volt
constexpr Real GEV = 1.0e9 * EV;      // Giga-Electron-Volt
constexpr Real JY = 1.e-23;           // Jansky
constexpr Real PC = 3.085678e18;      // Parsec
constexpr Real AU = 1.49597870691e13; // Astronomical unit
constexpr Real YEAR = 31536000.;
constexpr Real DAY = 86400.;
constexpr Real HOUR = 3600.;
constexpr Real MSUN = 1.989e33; // Solar mass
constexpr Real RSUN = 6.96e10;  // Solar radius
constexpr Real LSUN = 3.827e33; // Solar luminosity

} // namespace cgs
} // namespace constants
} // namespace singularity

#endif // SINGULARITY_OPAC_BASE_PHYSICAL_CONSTANTS_
