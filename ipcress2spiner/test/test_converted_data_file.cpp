//======================================================================
// ipcress2spiner tool for converting ipcresse to spiner
// Author: Alex Long (along@lanl.gov)
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

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <array>

#include <hdf5.h>
#include <hdf5_hl.h>

#ifndef SPINER_USE_HDF
#error "HDF5 must be enabled"
#endif

#include <ports-of-call/portability.hpp>
#include <singularity-opac/base/sp5.hpp>
#include <utils/spiner/spiner/sp5.hpp>
#include <utils/spiner/spiner/databox.hpp>
#include <utils/spiner/spiner/interpolation.hpp>

int main() {

  const char * filename = "test_file.h5";
  const char * group_name = "/10001";
  //const char * dataset_name = "/10001/plank total gray opacity";
  herr_t status = H5_SUCCESS;

  hid_t file;
  hid_t dataset;
  hid_t dataspace;

  // Open HDF5 file read-only
  file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  // Open dataset
  dataset = H5Gopen2(file, group_name, H5P_DEFAULT);

  // test gray databox
  {
    std::cout<<"--- Testing gray databox interpolation--- "<<std::endl;
    Spiner::DataBox<double> gray_databox;
    gray_databox.loadHDF(dataset, SP5::Fields::rgray);

    // make sure gray databox has the expected values
    if (gray_databox.rank() != 2) {
      std::cout<<"Rank check fails! Rank is: "<<gray_databox.rank()<<std::endl;
      return 1;
    }

    auto i_size = gray_databox.dim(2);
    auto j_size = gray_databox.dim(1);
    if (i_size != 30 || j_size != 30) {
      std::cout<<"size arguements do not match, N_1 = 3 and N_2 =3 but got: N_1 = "<<i_size;
      std::cout<<" and N_2 = "<<j_size<<std::endl;
      return 1;
    }

    // print values
    /*
    std::cout<<"i size: "<<i_size<<", j size: "<<j_size<<" total size: "<<gray_databox.size()<<std::endl;
    std::cout<<"element   value"<<std::endl;
    for (size_t j=0;j<j_size;++j) {
      for (size_t i=0;i<i_size;++i) {
        std::cout<<(j*i_size+i)<<"  "<<gray_databox(i,j)<<std::endl;
      }
    }
    */

    double interp_rho = log10(0.1); // g/cc

    double tolerance = 0.15;
    std::vector<double> T_points{0.1, 0.6, .99, 1.0, 1.01, 1.4, 3.0, 4.2, 6.2, 9.0, 10.0};
    // gold values from Draco's ipcress interpreter
    std::vector<double> gold_values{4.000100e-01, 4.758533e-01, 4.949213e-01, 0.5, 6.848793e-01,
                                    1.518317e+00, 1.879426e+01, 5.707131e+01, 2.064356e+02,
                                    7.064486e+02, 1.000307e+03};

    std::cout<<"T     gold    interpolated_value   relative_diff"<<std::endl;
    for (size_t i=0; i<T_points.size();++i) {
      double interp_T = log10(T_points[i]);
      double interped_value = gray_databox.interpToReal(interp_rho, interp_T);
      double relative_diff = std::fabs(interped_value-gold_values[i])/gold_values[i];
      std::cout<<std::pow(10,interp_T)<<"  "<<gold_values[i]<<"  "<<interped_value<<"  "<< relative_diff<<std::endl;
      if( relative_diff > tolerance) {
        std::cout<<"Value doesn't match to "<<tolerance<<"-- gold: "<<gold_values[i]<<" test: "<<interped_value<<std::endl;
        return 1;
      }
    }
  }

  // test multigroup databox
  {
    std::cout<<"--- Testing multigroup databox interpolation ---"<<std::endl;
    Spiner::DataBox<double> mg_databox;
    mg_databox.loadHDF(dataset, SP5::Fields::rsmg);

    // make sure mg databox has the expected values
    if (mg_databox.rank() != 3) {
      std::cout<<"Rank check fails! Rank is: "<<mg_databox.rank()<<std::endl;
      return 1;
    }

    auto i_size = mg_databox.dim(3);
    auto j_size = mg_databox.dim(2);
    auto k_size = mg_databox.dim(1);
    if (i_size != 12 || j_size != 6 ||  k_size != 6) {
      std::cout<<"size arguements do not match, N_1 = 3, N_2 = 3, N_3 = 3 but got: N_1 = "<<i_size;
      std::cout<<", N_2 = "<<j_size<<", and N_3 = "<<k_size<<std::endl;
      return 1;
    }
    std::cout<<"i size: "<<i_size<<", j size: "<<j_size<<" k_size: "<<k_size<<" total size: "<<mg_databox.size()<<std::endl;

    // print all values
    /*
    std::cout<<"element   value"<<std::endl;
    for (size_t k=0;k<k_size;++k) {
      for (size_t j=0;j<j_size;++j) {
        for (size_t i=0;i<i_size;++i) {
          std::cout<<(k*j_size*i_size + i_size*j + i)<<"  "<<mg_databox(k, j, i)<<std::endl;
        }
      }
    }
    */

    // interpolate at rho = 0.9 and T = 0.67 keV
    double interp_rho = log10(0.9); // g/cc
    double interp_T = log10(0.67); // keV

    double tolerance = 0.05;
    std::vector<double> hnu_points{0.02, 0.05, 0.085, 0.2, 0.5, 0.85, 2.0 , 5.0 , 8.5, 20.0, 50.0, 85.0};
    // gold values from Draco's ipcress interpreter
    std::vector<double> gold_values{0.3999584351001695, 0.39989004841500875, 0.399799302343056,
                                    0.39958473935770206, 0.39890319754138404, 0.39800204558786706,
                                    0.39588583373187775, 0.3892961262035159, 0.38087997690498276,
                                    0.3623440968644833, 0.3137366760288401, 0.2663122909516335};
    std::cout<<"hnu     gold    interpolated_value   relative_diff"<<std::endl;
    for (size_t i=0; i<hnu_points.size();++i) {
      double interp_hnu = log10(hnu_points[i]);
      double interped_value = mg_databox.interpToReal(interp_hnu, interp_rho, interp_T);
      double relative_diff = std::fabs(interped_value-gold_values[i])/gold_values[i];
      std::cout<<std::pow(10,interp_hnu)<<"  "<<gold_values[i]<<"  "<<interped_value<<"  "<< relative_diff<<std::endl;
      if( relative_diff > tolerance) {
        std::cout<<"Value doesn't match to "<<tolerance<<"-- gold: "<<gold_values[i]<<" test: "<<interped_value<<std::endl;
        return 1;
      }
    }
  }

  return (status == H5_SUCCESS) ? 0 : 1;
}
