// ======================================================================
// © 2026. Triad National Security, LLC. All rights reserved.  This
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

// This file was made in part with generative AI.

// Drives the SP5 table consistency checks by exit code: "good" must succeed,
// every other mode must fail.

#include <cstddef>
#include <iostream>
#include <string>

#include <hdf5.h>

#include <spiner/databox.hpp>

#include <singularity-opac/base/sp5.hpp>
#include <singularity-opac/photons/mean_opacity_photons.hpp>
#include <singularity-opac/photons/mean_s_opacity_photons.hpp>

using namespace singularity;

namespace {

using DataBox = Spiner::DataBox<Real>;

constexpr int NRho = 2;
constexpr int NT = 2;
constexpr int NGroups = 3;
constexpr char MaterialName[] = "guarded-material";

// Write a material group holding one rank-3 table of shape (nrho, nT, NGroups)
// and a group bounds dataset carrying nbounds entries. A consistent file has
// nrho == NRho, nT == NT, and nbounds == NGroups + 1. Bounds are strictly
// increasing and finite for every nbounds, so only the count is ever wrong.
void WriteTable(const std::string &filename, const char *field, const int nrho,
                const int nT, const int nbounds) {
  DataBox table(nrho, nT, NGroups);
  table.setRange(1, 2., 8., nT);
  table.setRange(2, -4., 2., nrho);
  for (int i = 0; i < table.size(); ++i) {
    table(i) = -2.;
  }

  DataBox bounds(nbounds);
  for (int i = 0; i < nbounds; ++i) {
    bounds(i) = 1.e11 * (i + 1);
  }

  hid_t file =
      H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  hid_t material =
      H5Gcreate(file, MaterialName, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  table.saveHDF(material, field);
  bounds.saveHDF(material, SP5::Multigroup::GroupBounds);
  H5Gclose(material);
  H5Fclose(file);

  table.finalize();
  bounds.finalize();
}

template <typename Mean>
int LoadAndCheck(const std::string &filename) {
  Mean loaded(filename, std::string(MaterialName));
  const int ngroups = loaded.ngroups();
  loaded.Finalize();
  if (ngroups != NGroups) {
    std::cerr << "unexpected ngroups: " << ngroups << std::endl;
    return 2;
  }
  return 0;
}

} // namespace

int main(int argc, char *argv[]) {
  const std::string usage =
      "usage: {absorption|scattering}-{good|short|long|flat_rho|flat_T}\n"
      "  good     bounds count matches the table, and must load\n"
      "  short    too few bounds; walking the groups would read past the end\n"
      "  long     too many bounds; every value read is valid, so only the "
      "count check can reject it\n"
      "  flat_rho a single density point; interpolation has no cell to "
      "bracket\n"
      "  flat_T   a single temperature point, likewise";
  if (argc != 2) {
    std::cerr << usage << std::endl;
    return 2;
  }
  const std::string mode = argv[1];
  const std::size_t split = mode.find('-');
  if (split == std::string::npos) {
    std::cerr << usage << std::endl;
    return 2;
  }
  const std::string reader = mode.substr(0, split);
  const std::string guard = mode.substr(split + 1);

  bool scattering;
  if (reader == "scattering") {
    scattering = true;
  } else if (reader == "absorption") {
    scattering = false;
  } else {
    std::cerr << usage << std::endl;
    return 2;
  }

  int nrho = NRho;
  int nT = NT;
  int nbounds = NGroups + 1;
  if (guard == "good") {
    // Nothing to perturb.
  } else if (guard == "short") {
    nbounds = NGroups - 1;
  } else if (guard == "long") {
    nbounds = NGroups + 3;
  } else if (guard == "flat_rho") {
    nrho = 1;
  } else if (guard == "flat_T") {
    nT = 1;
  } else {
    std::cerr << usage << std::endl;
    return 2;
  }

  const char *field = scattering ? SP5::MultigroupSOpac::RosselandGroupSOpacity
                                 : SP5::MultigroupOpac::RosselandGroupOpacity;
  const std::string filename = "sp5-table-guards-" + mode + ".sp5";
  WriteTable(filename, field, nrho, nT, nbounds);

  const int status = scattering
                         ? LoadAndCheck<photons::MeanSOpacityBase>(filename)
                         : LoadAndCheck<photons::MeanOpacityBase>(filename);
  if (status != 0) return status;

  std::cout << "loaded " << filename << std::endl;
  return 0;
}
