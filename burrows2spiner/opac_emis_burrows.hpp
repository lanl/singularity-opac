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

// Reader for Tabulated Neutrino Opacities                                    
// First published in Burrows, Reddy, Thompson, arXiv:astro-ph/0404432, 2004 
// Code from nubhlight code https://github.com/lanl/nubhlight       
// Miller, J. M., Ryan, B. R., Dolence, J. C. 2019, ApJS, 241:2
// Open-sourced under the same license

#ifndef BURROWS2SPINER_OPAC_EMIS_BURROWS_HPP_
#define BURROWS2SPINER_OPAC_EMIS_BURROWS_HPP_

#include <ports-of-call/portability.hpp>
#include <singularity-opac/base/radiation_types.hpp>

using singularity::RadiationType;
using singularity::RadType2Idx;
using singularity::Idx2RadType;

struct of_microphysics {
  double rho;
  double T;
  double Ye;
};

void   init_opac_emis_burrows();
void   fill_opac_emis_burrows(const Real rho, const Real T, const Real Ye);
double jnu_burrows(const double nu, const RadiationType type, const struct of_microphysics *m);
double Jnu_burrows(const double nu, const RadiationType type, const struct of_microphysics *m);
double int_jnudnudOmega_burrows(const struct of_microphysics *m);
double alpha_nu_burrows(double nu, int type, const struct of_microphysics *m);
void   test_opac_emis();

#endif // BURROWS2SPINER_OPAC_EMIS_BURROWS_HPP_
