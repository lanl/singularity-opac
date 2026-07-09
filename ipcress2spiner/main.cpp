//======================================================================
// ipcress2spiner tool for converting eospac to spiner
// Author: Alex Long (along@lanl.gov)
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

//#include "generate_files.hpp"
//#include "io_eospac.hpp"
#include "parse_cli.hpp"

/*
 *    Function prototypes for Gandolf C function equivalents
 */

extern "C" void c_gmatids( char*, int*,   int*,   int*,  int*);
extern "C" void     c_gkeys( char*, int*,char*[],   int*,  int*,  int*);
extern void   c_gsizeof( char*, int*,  char*,   int*,  int*);
extern "C" void  c_gchgrids( char*, int*,   int*,   int*,  int*,  int*,  int*, int*);
extern void  c_ggetgray( char*, int*,  char*,double*,  int*,  int*,
                       double*, int*,   int*,double*,  int*,  int*,  int*);
extern "C" void    c_ggetmg( char*, int*,  char*,double*,  int*,  int*,
                       double*, int*,   int*,
                       double*, int*,   int*,double*,  int*,  int*,  int*);
extern void  c_ggetdata( char*, int*,  char*,double*,  int*,  int*,  int*);
extern void  c_ggetchar( char*, int*,  char*,double*,  int*,  int*,  int*);

/*
 *    Function prototypes for Gandolf Fortran subroutines
 */

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

double interpolate_mg_opacity_data(const std::vector<double> &T_grid, const std::vector<double> &rho_grid, const std::vector<double> &group_bounds, const std::vector<double>  &op_data,
  double target_T, double target_rho, double target_hnu) {

  size_t n_T = T_grid.size();
  size_t n_rho = rho_grid.size();
  size_t n_groups = group_bounds.size()-1;

  auto T_L = std::distance(T_grid.begin(), std::lower_bound(T_grid.begin(), T_grid.end(), target_T));
  T_L = (target_T < T_grid[T_L]) ? T_L-1 : T_L;
  T_L = (T_L == (T_grid.size()-1)) ? T_L-1: T_L;
  auto T_R = T_L + 1;
  auto log_T_L = log(T_grid[T_L]);
  auto log_T_R = log(T_grid[T_R]);
  auto log_T_target = log(target_T);
  // fraction of T_2 to use
  auto T_frac = (log_T_target - log_T_L) / (log_T_R - log_T_L);

  auto rho_L = std::distance(rho_grid.begin(), std::lower_bound(rho_grid.begin(), rho_grid.end(), target_rho));
  rho_L = (target_rho < rho_grid[rho_L]) ? rho_L-1 : rho_L;
  rho_L = (rho_L == (rho_grid.size()-1)) ? rho_L-1: rho_L;
  auto rho_R = rho_L + 1;
  auto log_rho_L = log(rho_grid[rho_L]);
  auto log_rho_R = log(rho_grid[rho_R]);
  auto log_rho_target = log(target_rho);
  // fraction of rho_2 to use
  auto rho_frac = (log_rho_target - log_rho_L) / (log_rho_R - log_rho_L);

  auto group = std::distance(group_bounds.begin(), std::lower_bound(group_bounds.begin(), group_bounds.end(), target_hnu));
  group = (target_hnu < group_bounds[group]) ? group-1 : group;
  group = (group == (group_bounds.size()-1)) ? group-1 : group;

  //std::cout<<"Targets--T: "<<target_T<<" rho: "<<target_rho<<" hnu: "<<target_hnu<<std::endl;

  // get the adjacent rows of the opacity index, looks like op_[T]_[rho]_[hnu] with L indicating
  // left or index n and "R" indicating right index (n+1)
  auto op_L_L = op_data[  T_L * (n_rho*n_groups) + rho_L*(n_groups) + group];
  auto op_L_R = op_data[  T_L * (n_rho*n_groups) + rho_R*(n_groups) + group];
  auto op_R_L = op_data[  T_R * (n_rho*n_groups) + rho_L*(n_groups) + group];
  auto op_R_R = op_data[  T_R * (n_rho*n_groups) + rho_R*(n_groups) + group];

  // reduce in temperature dimension
  auto op_L =(1.0-T_frac) * op_L_L + T_frac*op_R_L;
  auto op_R =(1.0-T_frac) * op_L_R + T_frac*op_R_R;

  // reduce in density dimension and return
  double interp_opac = (1.0-rho_frac) * op_L + rho_frac*op_R;
  if (false) {
    std::cout<<"Interpolation box: "<<std::endl;
    std::cout<<op_L_R<<"--------"<<op_R_R<<" |"<<std::endl;
    std::cout<<"|          |"<<std::endl;
    std::cout<<"|          |"<<std::endl;
    std::cout<<"|  "<<interp_opac<<"     |"<<std::endl;
    std::cout<<"|          |"<<std::endl;
    std::cout<<"|          |"<<std::endl;
    std::cout<<op_L_L<<"--------"<<op_R_L<<std::endl;
  }
  return interp_opac;
}

Spiner::DataBox<double> build_opacity_spiner_databox(const std::vector<double> &temperature_points, const std::vector<double> &density_points, const std::vector<double> &group_bounds, const std::vector<double> &ramg) {

    Spiner::RegularGrid1D<double> new_temperature(temperature_points.front(), temperature_points.back(), temperature_points.size());
    Spiner::RegularGrid1D<double> new_density(density_points.front(), density_points.back(), density_points.size());
    Spiner::RegularGrid1D<double> new_group_bounds(group_bounds.front(), group_bounds.back(), 2*group_bounds.size());
    Spiner::DataBox<double> opacity_databox(temperature_points.size(), density_points.size(), 2*group_bounds.size());
    opacity_databox.setRange(0, new_temperature);
    opacity_databox.setRange(1, new_density);
    opacity_databox.setRange(2, new_group_bounds);

    std::cout<<"T     rho   hnu      opac(sq_cm/g)"<<std::endl;
    size_t flat_index=0;
    for (size_t k=0; k<new_temperature.nPoints();++k) {
      for (size_t j=0; j<new_density.nPoints();++j) {
        for (size_t i=0; i<(new_group_bounds.nPoints()-1);++i) {
          double new_group_center = 0.5*(new_group_bounds.x(i) + new_group_bounds.x(i+1));
          double interp_op = interpolate_mg_opacity_data(temperature_points, density_points, group_bounds, ramg, new_temperature.x(k), new_density.x(j), new_group_center);
          opacity_databox(k,j,i) = interp_op;
          std::cout<<new_temperature.x(k)<<"   "<<new_density.x(j)<<"    "<<new_group_center<<"     "<<opacity_databox(k,j,i)<<std::endl;
        }
      }
    }
  return opacity_databox;
}


int main(int argc, char *argv[]) {

  // ARL: Need these same parameters

  std::string filename;
  std::string savename, helpMessage;
  //Verbosity eospacWarn = Verbosity::Quiet;
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
  for (int imat=0;imat<nmats;++imat) {
    int mat_ID = mat_ids[imat];
    // return values from c_gchrids, initialized to invalid values
    int nt = -100;
    int nrho = -100;
    int nhnu = -100;
    int ngray = -100;
    int nmg = -100;
    // TODO ARL: Delete these...
    std::vector<char *> keys;
    for(int ikey = 0;ikey<32;++ikey) {
      keys.push_back(new char(16));
    }
    int kkeys = static_cast<int>(keys.size());
    int nkeys = -1;
    c_gkeys(filename_char, &mat_ID, keys.data(), &kkeys, &nkeys, &ier);
    std::cout<<"nkeys: "<<nkeys<<std::endl;
    std::vector<std::string> nice_keys;
    for (int j=0;j<nkeys;++j) {
      int key_length = 0;
      for (int k=0;k<16;++k) {
        if (keys[j][k] != '\0') key_length++;
      }
      nice_keys.push_back(std::string(keys[j], key_length));
    }
    for (auto ikey : nice_keys) {
      std::cout<<ikey<<" ";
    }
    std::cout<<std::endl;

    c_gchgrids(filename_char, &mat_ID, &nt, &nrho, &nhnu, &ngray, &nmg, &ier);

    print_gchrids_output(filename, mat_ID, nt, nrho, nhnu, ngray, nmg);

    std::vector<double> temperature_points(nt, 0.0);
    std::vector<double> group_bounds(nhnu, 0.0);
    std::vector<double> density_points(nrho, 0.0);
    int post_nt = nt;
    int post_nrho = nrho;
    int post_nhnu = nhnu;
    int post_nmg = nmg;

    std::string ramg_keyword = "ramg";
    std::string rsmg_keyword = "rsmg";
    std::string pmg_keyword = "pmg";

    std::array<std::string, 3> mg_opac_keywords{ramg_keyword, rsmg_keyword, pmg_keyword};

    for (auto key : mg_opac_keywords) {

      std::cout<<"Interpolating and building databox for "<<key<<std::endl;
      std::vector<double> opacity_data(nmg, 0.0);
      c_ggetmg(filename_char, &mat_ID, key.data(),
        temperature_points.data(), &nt, &post_nt,
        density_points.data(), &nrho, &post_nrho,
        group_bounds.data(), &nhnu, &post_nhnu,
        opacity_data.data(), &nmg, &post_nmg, &ier);

      auto spiner_opacity_databox = build_opacity_spiner_databox(temperature_points, density_points, group_bounds, opacity_data);
      // save here with keyword?
    }
  }

  std::cout << "Done." << std::endl;

  return (status == H5_SUCCESS) ? 0 : 1;
}
