//======================================================================
// spiner_databox.h make spiner databox structure from ipcress data
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

#ifndef _IPCRESS2SPINER_MAKE_SPINER_DATABOX_HPP_
#define _IPCRESS2SPINER_MAKE_SPINER_DATABOX_HPP_

inline std::pair<Spiner::DataBox<double>, std::vector<double>>  build_multigroup_opacity_spiner_databox(const std::vector<double> &temperature_points, const std::vector<double> &density_points, const std::vector<double> &group_bounds, const std::vector<double> &full_multigroup_data, const bool log_T_rho_hnu, const bool log_interpolation, const bool verbose_mode = false) {

  const auto n_groups = group_bounds.size() -1;
  const double temperature_start = (log_T_rho_hnu) ? log10(temperature_points.front()) : temperature_points.front();
  const double temperature_end = (log_T_rho_hnu) ? log10(temperature_points.back()) : temperature_points.back();
  const double density_start = (log_T_rho_hnu) ? log10(density_points.front()) : density_points.front();
  const double density_end = (log_T_rho_hnu) ? log10(density_points.back()) : density_points.back();
  const double hnu_start = (log_T_rho_hnu) ? log10(group_bounds.front()) : group_bounds.front();
  const double hnu_end = (log_T_rho_hnu) ? log10(group_bounds.back()) : group_bounds.back();
  const double hnu_midpoints_start = (log_T_rho_hnu) ? log10( 0.5*(group_bounds[0]+group_bounds[1])) : 0.5*(group_bounds[0] + group_bounds[1]);
  const double hnu_midpoints_end = (log_T_rho_hnu) ? log10(0.5*(group_bounds[n_groups]+group_bounds[n_groups-1])) : 0.5*(group_bounds[n_groups]+group_bounds[n_groups-1]);

  Spiner::RegularGrid1D<double> new_temperature(temperature_start, temperature_end, temperature_points.size());
  Spiner::RegularGrid1D<double> new_density(density_start, density_end, density_points.size());
  Spiner::RegularGrid1D<double> new_group_bounds(hnu_start, hnu_end, group_bounds.size());
  Spiner::RegularGrid1D<double> new_group_centers(hnu_midpoints_start, hnu_midpoints_end, group_bounds.size()-1);
  Spiner::DataBox<double> opacity_databox(temperature_points.size(), density_points.size(), n_groups);
  opacity_databox.setRange(0, new_temperature);
  opacity_databox.setRange(1, new_density);
  opacity_databox.setRange(2, new_group_centers);

  if(verbose_mode) {
    std::cout<<"T     rho   hnu      opac(sq_cm/g)"<<std::endl;
  }

  size_t flat_index=0;
  for (size_t k=0; k<new_temperature.nPoints();++k) {
    for (size_t j=0; j<new_density.nPoints();++j) {
      for (size_t i=0; i<(new_group_bounds.nPoints()-1);++i) {
        double new_group_center = 0.5*(new_group_bounds.x(i) + new_group_bounds.x(i+1));
        if (log_T_rho_hnu)
          new_group_center = 0.5*(std::pow(10, new_group_bounds.x(i)) + std::pow(10, new_group_bounds.x(i+1)));
        double target_T = (log_T_rho_hnu) ? std::pow(10, new_temperature.x(k)) : new_temperature.x(k);
        double target_rho = (log_T_rho_hnu) ? std::pow(10, new_density.x(j)) : new_density.x(j);
        double interp_op = interpolate_mg_opacity_data(temperature_points, density_points, group_bounds, full_multigroup_data, target_T, target_rho, new_group_center, log_interpolation);
        opacity_databox(k,j,i) = interp_op;
        if(verbose_mode) {
          std::cout<<target_T <<"   "<<target_rho<<"    "<<new_group_center<<"     "<<opacity_databox(k,j,i)<<std::endl;
        }
      }
    }
  }

  // return new group bounds to write separately to hdf5 file
  std::vector<double> vector_new_group_bounds(group_bounds.size());
  for (size_t i=0;i<group_bounds.size();++i) {
    vector_new_group_bounds[i]= new_group_bounds.x(i);
  }

  return {opacity_databox, vector_new_group_bounds};
}

inline Spiner::DataBox<double> build_gray_opacity_spiner_databox(const std::vector<double> &temperature_points, const std::vector<double> &density_points, const std::vector<double> &full_gray_data, const bool log_T_rho_hnu, const bool log_interpolation, const bool verbose_mode = false) {

  const double temperature_start = (log_T_rho_hnu) ? log10(temperature_points.front()) : temperature_points.front();
  const double temperature_end = (log_T_rho_hnu) ? log10(temperature_points.back()) : temperature_points.back();
  const double density_start = (log_T_rho_hnu) ? log10(density_points.front()) : density_points.front();
  const double density_end = (log_T_rho_hnu) ? log10(density_points.back()) : density_points.back();

  Spiner::RegularGrid1D<double> new_temperature(temperature_start, temperature_end, temperature_points.size());
  Spiner::RegularGrid1D<double> new_density(density_start, density_end, density_points.size());
  Spiner::DataBox<double> opacity_databox(temperature_points.size(), density_points.size());
  opacity_databox.setRange(0, new_temperature);
  opacity_databox.setRange(1, new_density);

  if(verbose_mode) {
    std::cout<<"T     rho   opac(sq_cm/g)"<<std::endl;
  }

  size_t flat_index=0;
  for (size_t k=0; k<new_temperature.nPoints();++k) {
    for (size_t j=0; j<new_density.nPoints();++j) {
      double target_T = (log_T_rho_hnu) ? std::pow(10, new_temperature.x(k)) : new_temperature.x(k);
      double target_rho = (log_T_rho_hnu) ? std::pow(10, new_density.x(j)) : new_density.x(j);
      double interp_op = interpolate_gray_opacity_data(temperature_points, density_points, full_gray_data, target_T, target_rho, log_interpolation);
      opacity_databox(k,j) = interp_op;
      if(verbose_mode) {
        std::cout<<target_T <<"   "<<target_rho<<"    "<<opacity_databox(k,j)<<std::endl;
      }
    }
  }

  return opacity_databox;
}

#endif // _IPCRESS2SPINER_MAKE_SPINER_DATABOX_HPP_

