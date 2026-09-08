//======================================================================
// generate_files function that write spiner file
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

#ifndef _IPCRESS2SPINER_GENERATE_FILES_HPP_
#define _IPCRESS2SPINER_GENERATE_FILES_HPP_

#include <string>
#include <vector>
#include <utils/spiner/spiner/sp5.hpp>
#include <utils/spiner/spiner/databox.hpp>
#include <utils/spiner/spiner/interpolation.hpp>

herr_t saveMaterial(hid_t loc, hid_t matGroup, const int matid,
                    const std::string &sMatid,
                    const std::string &sp5_field_name,
                    Spiner::DataBox<double> &opacity) {


  double zero_offset =0.0;
  herr_t status = 0;
  // Dependent variables metadata
  status += H5LTset_attribute_string(loc, sMatid.c_str(), SP5::Offsets::messageName,
                                     SP5::Offsets::message);
  status += H5LTset_attribute_double(loc, sMatid.c_str(), SP5::Offsets::rho,
                                     &zero_offset, 1);
  status += H5LTset_attribute_double(loc, sMatid.c_str(), SP5::Offsets::T, &zero_offset, 1);

  status += opacity.saveHDF(matGroup, sp5_field_name);

  return status;
}

#endif // _IPCRESS2SPINER_GENERATE_FILES_HPP_
