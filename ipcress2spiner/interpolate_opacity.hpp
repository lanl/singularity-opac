//======================================================================
// ipcress2spiner tool for converting ipcress to spiner
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

#ifndef _IPCRESS2SPINER_INTERPOLATE_HPP_
#define _IPCRESS2SPINER_INTERPOLATE_HPP_

inline double interpolate_mg_opacity_data(const std::vector<double> &T_grid, const std::vector<double> &rho_grid, const std::vector<double> &group_bounds, const std::vector<double>  &op_data,
  double target_T, double target_rho, double target_hnu) {

  size_t n_T = T_grid.size();
  size_t n_rho = rho_grid.size();
  size_t n_groups = group_bounds.size()-1;

  // always use log interpolation between temperatures, densities, and opacity table values
  auto T_L = std::distance(T_grid.begin(), std::lower_bound(T_grid.begin(), T_grid.end(), target_T));
  T_L = (target_T < T_grid[T_L]) ? T_L-1 : T_L;
  T_L = (T_L == (T_grid.size()-1)) ? T_L-1: T_L;
  auto T_R = T_L + 1;
  // fraction of T_R opacity to use
  double log_T_L = log(T_grid[T_L]);
  double log_T_R = log(T_grid[T_R]);
  double log_T_target = log(target_T);
  double T_frac = (log_T_target - log_T_L) / (log_T_R - log_T_L);

  auto rho_L = std::distance(rho_grid.begin(), std::lower_bound(rho_grid.begin(), rho_grid.end(), target_rho));
  rho_L = (target_rho < rho_grid[rho_L]) ? rho_L-1 : rho_L;
  rho_L = (rho_L == (rho_grid.size()-1)) ? rho_L-1: rho_L;
  auto rho_R = rho_L + 1;
  // fraction of rho_R opacity to use
  double log_rho_L = log(rho_grid[rho_L]);
  double log_rho_R = log(rho_grid[rho_R]);
  double log_rho_target = log(target_rho);
  double rho_frac = (log_rho_target - log_rho_L) / (log_rho_R - log_rho_L);

  auto group = std::distance(group_bounds.begin(), std::lower_bound(group_bounds.begin(), group_bounds.end(), target_hnu));
  group = (target_hnu < group_bounds[group]) ? group-1 : group;
  group = (group == (group_bounds.size()-1)) ? group-1 : group;

  // get the adjacent rows of the opacity index, looks like op_[T]_[rho]_[hnu] with L indicating
  // left or index n and "R" indicating right index (n+1)
  double op_L_L = log(op_data[  T_L * (n_rho*n_groups) + rho_L*(n_groups) + group]);
  double op_L_R = log(op_data[  T_L * (n_rho*n_groups) + rho_R*(n_groups) + group]);
  double op_R_L = log(op_data[  T_R * (n_rho*n_groups) + rho_L*(n_groups) + group]);
  double op_R_R = log(op_data[  T_R * (n_rho*n_groups) + rho_R*(n_groups) + group]);

  // reduce in temperature dimension
  auto op_L =op_L_L + T_frac*(op_R_L - op_L_L);
  auto op_R =op_L_R + T_frac*(op_R_R - op_L_R);

  // reduce in density dimension and return
  double interp_opac = exp(op_L + rho_frac*(op_R -op_L));
  if (false) {
    std::cout<<"Targets--T: "<<target_T<<" rho: "<<target_rho<<" hnu: "<<target_hnu<<" T frac:" <<T_frac<<" rho frac: "<<rho_frac<<std::endl;
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

inline double interpolate_gray_opacity_data(const std::vector<double> &T_grid, const std::vector<double> &rho_grid, const std::vector<double>  &op_data,
  double target_T, double target_rho) {

  size_t n_T = T_grid.size();
  size_t n_rho = rho_grid.size();

  // always use log interpolation between temperatures, densities, and opacity table values
  auto T_L = std::distance(T_grid.begin(), std::lower_bound(T_grid.begin(), T_grid.end(), target_T));
  T_L = (target_T < T_grid[T_L]) ? T_L-1 : T_L;
  T_L = (T_L == (T_grid.size()-1)) ? T_L-1: T_L;
  auto T_R = T_L + 1;
  // fraction of T_R opacity to use
  double log_T_L = log(T_grid[T_L]);
  double log_T_R = log(T_grid[T_R]);
  double log_T_target = log(target_T);
  double T_frac = (log_T_target - log_T_L) / (log_T_R - log_T_L);

  auto rho_L = std::distance(rho_grid.begin(), std::lower_bound(rho_grid.begin(), rho_grid.end(), target_rho));
  rho_L = (target_rho < rho_grid[rho_L]) ? rho_L-1 : rho_L;
  rho_L = (rho_L == (rho_grid.size()-1)) ? rho_L-1: rho_L;
  auto rho_R = rho_L + 1;
  // fraction of rho_R opacity to use
  double log_rho_L = log(rho_grid[rho_L]);
  double log_rho_R = log(rho_grid[rho_R]);
  double log_rho_target = log(target_rho);
  double rho_frac = (log_rho_target - log_rho_L) / (log_rho_R - log_rho_L);

  // get the adjacent rows of the opacity index, looks like op_[T]_[rho]_[hnu] with L indicating
  // left or index n and "R" indicating right index (n+1)
  double op_L_L = log(op_data[  T_L * n_rho + rho_L]);
  double op_L_R = log(op_data[  T_L * n_rho + rho_R]);
  double op_R_L = log(op_data[  T_R * n_rho + rho_L]);
  double op_R_R = log(op_data[  T_R * n_rho + rho_R]);

  // reduce in temperature dimension
  double op_L =op_L_L + T_frac*(op_R_L - op_L_L);
  double op_R =op_L_R + T_frac*(op_R_R - op_L_R);

  // reduce in density dimension and return
  double interp_opac = exp(op_L + rho_frac*(op_R - op_L));
  if (false) {
    std::cout<<"Targets--T: "<<target_T<<" rho: "<<target_rho<< "T frac:" <<T_frac<<" rho frac: "<<rho_frac<<std::endl;
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

#endif // _IPCRESS2SPINER_INTERPOLATE_HPP_
