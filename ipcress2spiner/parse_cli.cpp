//======================================================================
// sesame2spiner tool for converting eospac to spiner
// Author: Jonah Miller (jonahm@lanl.gov)
// © 2021-2023. Triad National Security, LLC. All rights reserved.  This
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

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "parse_cli.hpp"

void parseCLI(int argc, char *argv[], std::string &savename,
              std::string &filename, bool &printMetadata,
              std::string &helpMessage) {

  std::stringstream helpStream;
  helpStream << "Usage: " << argv[0]
             << "[-p] [-w] [-h] [-v] [-vv] [-d] [-s <savename>] <parameter file>\n\n"
             << "\t <ipcress file>: input ipcress file\n"
             << "\t-s <savename>: filename to save to. Defaults to " << DEFAULT_SAVENAME
             << "\n"
             << "\t-p:  print metadata associated with materials "
             << "in parameter files\n"
             << "\t-v:  print ipcress warnings\n"
             << "\t-vv: print debug information\n"
             << "\t-w:  same as -v\n"
             << "\t-d:  same as -vv\n"
             << "\t-h:  print this message\n"
             << "\n"
             << std::endl;
  helpMessage = helpStream.str();

  savename = DEFAULT_SAVENAME;

  if (argc < 2) {
    std::cerr << helpMessage << std::endl;
    std::exit(1);
  }
  if (argc == 2 && std::strcmp(argv[1], "-h") == 0) {
    std::cout << helpMessage << std::endl;
    std::exit(0);
  }
  for (int i = 1; i < argc; i++) {
    if (std::strcmp(argv[i], "-h") == 0) {
      std::cout << helpMessage << std::endl;
      std::exit(0);
    } else if (std::strcmp(argv[i], "-p") == 0) {
      printMetadata = true;
    } else if (std::strcmp(argv[i], "-s") == 0) {
      savename = argv[++i];
    } else {
      filename = std::string(argv[i]);
    }
  }
}
