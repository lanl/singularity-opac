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
#include <gandolf.h>

#include "generate_files.hpp"
#include "interpolate_opacity.hpp"
#include "parse_cli.hpp"
#include "make_spiner_databox.hpp"

// ARL: My to-do list
// \TODO Thread equal log spacing of rho and T points through command line
// \TODO Thread verbosity option through command-line options to print more info
// \TODO Add ionization state? Leave until singularity-eos and singularity-opac are merged
// \TODO clang-format
// \TODO Delete char fields from new[]

// Function prototypes for Gandolf C function equivalents, prepend "C" to prevent name mangling
// gandolf.h does not expose the interface to these functions, which is why they are declared here
extern "C" void c_gmatids( char*, int*,   int*,   int*,  int*);
extern "C" void     c_gkeys( char*, int*,char*[],   int*,  int*,  int*);
extern "C" void   c_gsizeof( char*, int*,  char*,   int*,  int*);
extern "C" void  c_gchgrids( char*, int*,   int*,   int*,  int*,  int*,  int*, int*);
extern "C" void  c_ggetgray( char*, int*,  char*,double*,  int*,  int*,
                       double*, int*,   int*,double*,  int*,  int*,  int*);
extern "C" void    c_ggetmg( char*, int*,  char*,double*,  int*,  int*,
                       double*, int*,   int*,
                       double*, int*,   int*,double*,  int*,  int*,  int*);
extern void  c_ggetdata( char*, int*,  char*,double*,  int*,  int*,  int*);
extern void  c_ggetchar( char*, int*,  char*,double*,  int*,  int*,  int*);

// Function prototypes for Gandolf Fortran subroutines
// gandolf.h does not expose the interface to these functions, which is why they are declared here
extern void   f_gmatids( char*, int*,   int*,  int*,  int*);
extern void     f_gkeys( char*, int*,  char*,  int*,  int*,   int*);
extern void   f_gsizeof( char*, int*,  char*,  int*,  int*);
extern void  f_gchgrids( char*, int*,   int*,  int*,  int*,   int*,   int*, int*);
extern void  f_ggetgray( char*, int*,  char*,double*, int*,   int*,double*,
                          int*, int*,double*,   int*, int*,   int*);
extern void    f_ggetmg( char*, int*,  char*,double*, int*,   int*,double*,
                          int*, int*,double*,   int*, int*,double*,   int*, int*, int*);
extern void  f_ggetdata( char*, int*,  char*,double*, int*,   int*,   int*);
extern void  f_ggetchar( char*, int*,  char*,double*, int*,   int*,   int*);

//constexpr int H5_SUCCESS = 0; // older HDF5's declare this, just declare it here

void print_gchrids_output(const std::string &filename, const int mat_ID, const int nt, const int nrho, const int nhnu, const int ngray, const int nmg) {
  std::cout<<"---- Gandolf information for material "<<mat_ID<<" in "<<filename<<" ----"<<std::endl;
  std::cout<<"Temperature poitns: "<<nt<<" density points: "<<nrho<<" frequency points: "<<nhnu<<std::endl;
  std::cout<<"number of gray points: "<<ngray<<" number of multigroup points: "<<nmg<<std::endl;
}

void print_ggetmg_output(const int nmg, const int post_ramg_nmg,
    const std::vector<double> &temperature_points,
    const std::vector<double> &density_points,
    const std::vector<double> &group_bounds,
    const std::vector<double> &ramg) {

    std::cout<<"Temperature points: "<<temperature_points.size()<<std::endl;
    for (auto it : temperature_points)
      std::cout<<it<<" ";
    std::cout<<std::endl;

    std::cout<<"Density points: "<<density_points.size()<<std::endl;
    for (auto it : density_points)
      std::cout<<it<<" ";
    std::cout<<std::endl;

    std::cout<<"group bounds: "<<group_bounds.size()<<std::endl;
    for (auto it : group_bounds)
      std::cout<<it<<" ";
    std::cout<<std::endl;

    std::cout<<" Raw data, pre size: "<<nmg<<" post size: "<<post_ramg_nmg<<std::endl;
    for (size_t idata =0; idata<static_cast<size_t>(nmg); ++idata) {
      std::cout<<idata<<" "<<ramg[idata]<<std::endl;
    }
}

int main(int argc, char *argv[]) {

  std::string filename;
  std::string savename, helpMessage;
  bool printMetadata = false;
  herr_t status = H5_SUCCESS;

  parseCLI(argc, argv, savename, filename, printMetadata, helpMessage);

  std::cout << "ipcress2spiner                            \n"
            << "-----------------------------------------\n"
            << "Author: Alex Long (along@lanl.gov)   \n"
            << "-----------------------------------------\n"
            << std::endl;

  char * filename_char = filename.data();
  int ier = -100;

  std::cout<<"Getting material IDs from file: "<<filename<<std::endl;
  std::vector<int> mat_ids(256, -1);
  int kmats = static_cast<int>(mat_ids.size());
  int nmats = -1;
  c_gmatids(filename_char, mat_ids.data(), &kmats, &nmats, &ier);

  std::cout<<"nmats: "<<nmats<<std::endl;

  // output HDF5 setup
  std::cout << "Saving to file " << savename << std::endl;
  auto file_loc = H5Fcreate(savename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  H5LTset_attribute_string(file_loc, "/", "ipcress2spiner version", IPCRESS2SPINER_VERSION);

  for (int imat=0;imat<nmats;++imat) {
    int mat_ID = mat_ids[imat];
    std::cout<<"Processing material "<<mat_ID<<std::endl;
    // create this mat ID in the HDF5 file
    herr_t status = 0;
    auto sMatID = std::to_string(mat_ID);

    hid_t matGroup = H5Gcreate(file_loc, sMatID.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    // ARL: Is this softlink necessary?
    //status += H5Lcreate_soft(sMatID.c_str(), file_loc, savename.c_str(), H5P_DEFAULT, H5P_DEFAULT);
    if (status != H5_SUCCESS) {
      std::cerr << "WARNING: problem with HDf5 post H5GCreate call" << std::endl;
    }

    // return values from c_gchrids, initialized to invalid values
    int nt = -100;
    int nrho = -100;
    int nhnu = -100;
    int ngray = -100;
    int nmg = -100;
    // 32 keys has been enough for all keys in the files so far
    constexpr size_t total_keys = 32;
    constexpr size_t max_key_length = 16;
    std::vector<char *> keys;
    for (size_t ikey =0;ikey<total_keys;++ikey) {
      keys.push_back(new char[max_key_length]);
    }

    int kkeys = static_cast<int>(keys.size());
    int nkeys = -1;
    c_gkeys(filename_char, &mat_ID, keys.data(), &kkeys, &nkeys, &ier);
    if(nkeys > total_keys) {
      std::cout<<"POTENTIAL ERROR: nkeys: "<<nkeys<<" is larger than current max keys"<<total_keys<<std::endl;
    }
    std::vector<std::string> nice_keys;
    for (int j=0;j<nkeys;++j) {
      int key_length = 0;
      for (int k=0;k<max_key_length;++k) {
        if (keys[j][k] != '\0') key_length++;
      }
      nice_keys.push_back(std::string(keys[j], key_length));
    }

    // print the nice string form of the keys
    for (auto ikey : nice_keys) {
      std::cout<<ikey<<" ";
    }
    std::cout<<std::endl;

    c_gchgrids(filename_char, &mat_ID, &nt, &nrho, &nhnu, &ngray, &nmg, &ier);

    print_gchrids_output(filename, mat_ID, nt, nrho, nhnu, ngray, nmg);

    std::vector<double> temperature_points(nt, 0.0);
    std::vector<double> group_bounds(nhnu, 0.0);
    std::vector<double> density_points(nrho, 0.0);
    // these post values are potentially reset by calls to gandolf?
    int post_nt = nt;
    int post_nrho = nrho;
    int post_nhnu = nhnu;
    int post_nmg = nmg;
    int post_ngray = ngray;

    std::string ramg_keyword = "ramg";
    std::string rsmg_keyword = "rsmg";
    std::string rtmg_keyword = "rtmg";
    std::string pmg_keyword = "pmg";
    std::string ragray_keyword = "ragray";
    std::string rgray_keyword = "rgray";
    std::string pgray_keyword = "pgray";

    std::array<std::string, 7> mg_opac_keywords{ramg_keyword, rsmg_keyword, rtmg_keyword, pmg_keyword,
                                                ragray_keyword, rgray_keyword, pgray_keyword};
    std::array<std::string, 7> mg_fields{ SP5::Fields::ramg, SP5::Fields::rsmg, SP5::Fields::rtmg, SP5::Fields::pmg,
                                          SP5::Fields::ragray, SP5::Fields::rgray, SP5::Fields::pgray};

    // TODO Thread these variables through the input
    bool log_T_rho_hnu = true; // form a new grid from max and min evenly spaced in log10

    for (size_t i=0;i<mg_opac_keywords.size();++i) {
      auto key = mg_opac_keywords[i];
      auto sp5_field_name = mg_fields[i];
      std::cout<<"Interpolating and building databox for "<<mat_ID<<" and "<<key<<std::endl;
      // multigroup opacities
      if(key == "ramg" || key == "rsmg" || key=="rtmg" || key == "pmg") {
        std::vector<double> opacity_data(nmg, 0.0);
        c_ggetmg(filename_char, &mat_ID, key.data(),
          temperature_points.data(), &nt, &post_nt,
          density_points.data(), &nrho, &post_nrho,
          group_bounds.data(), &nhnu, &post_nhnu,
          opacity_data.data(), &nmg, &post_nmg, &ier);

        auto [spiner_opacity_databox, new_group_bounds] = build_multigroup_opacity_spiner_databox(
          temperature_points, density_points, group_bounds, opacity_data, log_T_rho_hnu);

        // only save the group bounds once for this material
        if (i==0) {
          const hsize_t dimensions[] = {static_cast<hsize_t>(new_group_bounds.size())};
          std::string dataset_name("group bounds");
          hid_t dataspace_id = H5Screate_simple(1, dimensions, nullptr);
          hid_t dataset_id = H5Dcreate2(matGroup, dataset_name.c_str(), H5T_IEEE_F64LE, dataspace_id,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
          const herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
              new_group_bounds.data());
        }
        std::cout<<"Saving multigroup databox for "<<mat_ID<<" and "<<key<<std::endl;
        saveMaterial(file_loc, matGroup, mat_ID, sMatID, sp5_field_name, spiner_opacity_databox);
      }
      // gray opacites
      else {
        std::vector<double> opacity_data(ngray, 0.0);
       c_ggetgray(filename_char, &mat_ID, key.data(),
          temperature_points.data(), &nt, &post_nt,
          density_points.data(), &nrho, &post_nrho,
          opacity_data.data(), &ngray, &post_ngray, &ier);

        auto spiner_opacity_databox = build_gray_opacity_spiner_databox(
          temperature_points, density_points, opacity_data, log_T_rho_hnu);

        std::cout<<"Saving gray databox for "<<mat_ID<<" and "<<key<<std::endl;
        saveMaterial(file_loc, matGroup, mat_ID, sMatID, sp5_field_name, spiner_opacity_databox);
      }
    }
    status += H5Gclose(matGroup);
  }

  status += H5Fclose(file_loc);
  if (status != H5_SUCCESS) {
    std::cerr << "WARNING: problem with HDf5" << std::endl;
  }

  std::cout << "Done." << std::endl;

  return (status == H5_SUCCESS) ? 0 : 1;
}
