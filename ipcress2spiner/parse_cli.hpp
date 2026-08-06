//======================================================================
// ipcress2spiner tool for converting eospac to spiner
// Author: Alex R. Long (along@lanl.gov)
// © 2026. Triad National Security, LLC. All rights reserved.  This
// program was produced under U.S. Government contract 89233218CNA000001
// for Los Alamos National Laboratory (LANL), which is operated by Triad
// National Security, LLC for the U.S.  Department of Energy/National
// Nuclear Security Administration. All rights in the program are
// reserved by Triad National Security, LLC, and the U.S. Department of
// Energy/National Nuclear Security Administration. The Government is
// granted for itself and others acting on its behalf a nonexclusive,
// paid-up, irrevocable worldwide license in this material to reproduce,
// prepare derivative works, distribute copies to the public, perform
// publicly and display publicly, and to permit others to do so.
//======================================================================

#ifndef _IPCRESS2SPINER_PARSER_HPP_
#define _IPCRESS2SPINER_PARSER_HPP_

#include <string>
#include <vector>

const std::string DEFAULT_SAVENAME = "converted_ipcress.h5";;

void parseCLI(int argc, char *argv[], std::string &savename,
              std::string &filename, bool &printMetadata,
              std::string &helpMessage);

#endif // _IPCRESS2SPINER_PARSER_HPP_
