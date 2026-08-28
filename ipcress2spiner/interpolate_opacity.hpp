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


namespace {

  std::pair<size_t, size_t> get_bounding_indices(const std::vector<double> &grid, const double target_value) {
    auto left_index = std::distance(grid.begin(), std::lower_bound(grid.begin(), grid.end(), target_value));
    left_index = (target_value < grid[left_index]) ? left_index-1 : left_index;
    left_index = (left_index == (grid.size()-1)) ? left_index-1: left_index;
    auto right_index = left_index + 1;
    return {left_index, right_index};
  }

  double get_fraction(const std::vector<double> &grid, const size_t left_index, const size_t right_index, const double target_value) {
    // always use log interpolation between temperatures, densities, and opacity table values
    double log_left_value = log(grid[left_index]);
    double log_right_value = log(grid[right_index]);
    double log_target_value = log(target_value);
    // fraction is defined in terms of amount of right value to use in interpolation
    double fraction_right_value = (log_target_value - log_left_value) / (log_right_value - log_left_value);
    return fraction_right_value;
  }

  double interpolate_opacity_data(const size_t T_L, const size_t T_R, const size_t rho_L,
                          const size_t rho_R, const size_t group, const double T_frac,
                          const double rho_frac, const size_t n_rho, const size_t n_groups,
                          const std::vector<double> &op_data) {
    // get the adjacent rows of the opacity index, looks like op_[T]_[rho]_[hnu] with L indicating
    // left or index n and "R" indicating right index (n+1)
    double op_L_L = log(op_data[  T_L * (n_rho*n_groups) + rho_L*(n_groups) + group]);
    double op_L_R = log(op_data[  T_L * (n_rho*n_groups) + rho_R*(n_groups) + group]);
    double op_R_L = log(op_data[  T_R * (n_rho*n_groups) + rho_L*(n_groups) + group]);
    double op_R_R = log(op_data[  T_R * (n_rho*n_groups) + rho_R*(n_groups) + group]);

    // reduce in temperature dimension
    auto op_L =op_L_L + T_frac*(op_R_L - op_L_L);
    auto op_R =op_L_R + T_frac*(op_R_R - op_L_R);

    // reduce in density dimension
    double interp_opac = exp(op_L + rho_frac*(op_R -op_L));
    // get the adjacent rows of the opacity index, looks like op_[T]_[rho]_[hnu] with L indicating
    if (false) {
      //std::cout<<"Targets--T: "<<target_T<<" rho: "<<target_rho<<" hnu: "<<target_hnu<<" T frac:" <<T_frac<<" rho frac: "<<rho_frac<<std::endl;
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
}

inline double interpolate_mg_opacity_data(const std::vector<double> &T_grid, const std::vector<double> &rho_grid, const std::vector<double> &group_bounds, const std::vector<double>  &op_data,
  double target_T, double target_rho, double target_hnu) {

  size_t n_T = T_grid.size();
  size_t n_rho = rho_grid.size();
  size_t n_groups = group_bounds.size()-1;

  auto [T_L, T_R] = get_bounding_indices(T_grid, target_T);
  auto T_frac = get_fraction(T_grid, T_L, T_R, target_T);
  auto [rho_L, rho_R] = get_bounding_indices(rho_grid, target_rho);
  auto rho_frac = get_fraction(rho_grid, rho_L, rho_R, target_rho);

  auto group = std::distance(group_bounds.begin(), std::lower_bound(group_bounds.begin(), group_bounds.end(), target_hnu));
  group = (target_hnu < group_bounds[group]) ? group-1 : group;
  group = (group == (group_bounds.size()-1)) ? group-1 : group;

  return interpolate_opacity_data(T_L, T_R, rho_L, rho_R, group, T_frac, rho_frac, n_rho, n_groups, op_data);
}

inline double interpolate_gray_opacity_data(const std::vector<double> &T_grid, const std::vector<double> &rho_grid, const std::vector<double>  &op_data,
  double target_T, double target_rho) {

  size_t n_T = T_grid.size();
  size_t n_rho = rho_grid.size();

  auto [T_L, T_R] = get_bounding_indices(T_grid, target_T);
  auto T_frac = get_fraction(T_grid, T_L, T_R, target_T);
  auto [rho_L, rho_R] = get_bounding_indices(rho_grid, target_rho);
  auto rho_frac = get_fraction(rho_grid, rho_L, rho_R, target_rho);

  constexpr size_t group = 0; // always 0 for gray
  constexpr size_t n_groups = 1; // always 1 for gray
  return interpolate_opacity_data(T_L, T_R, rho_L, rho_R, group, T_frac, rho_frac, n_rho, n_groups, op_data);
}

#endif // _IPCRESS2SPINER_INTERPOLATE_HPP_
