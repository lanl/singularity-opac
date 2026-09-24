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
#ifndef SINGULARITY_OPAC_PHOTONS_MEAN_PHOTON_UTILS_
#define SINGULARITY_OPAC_PHOTONS_MEAN_PHOTON_UTILS_
// This file was made in part with generative AI.

#include <cmath>
#include <limits>

#include <ports-of-call/portability.hpp>
#include <singularity-opac/base/opac_error.hpp>
#include <singularity-opac/base/sp5.hpp>
#include <singularity-opac/photons/mean_photon_types.hpp>
#include <singularity-opac/photons/thermal_distributions_photons.hpp>
#include <spiner/databox.hpp>

#ifdef SPINER_USE_HDF
#include <filesystem>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <optional>
#include <string>
#endif

namespace singularity {
namespace photons {
namespace impl {

using MeanUtilsDataBox = Spiner::DataBox<Real>;

#ifdef SPINER_USE_HDF
// Sentinels returned by the material-group openers below.
constexpr hid_t MaterialNotFound = -1;
constexpr hid_t MaterialAmbiguous = -2;

// Suppresses HDF5's automatic error printing for the current scope and restores
// whatever handler was installed on the way out.
class ScopedH5ErrorHandler {
 public:
  ScopedH5ErrorHandler() {
    H5Eget_auto2(H5E_DEFAULT, &func_, &data_);
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
  }
  ~ScopedH5ErrorHandler() { H5Eset_auto2(H5E_DEFAULT, func_, data_); }
  ScopedH5ErrorHandler(const ScopedH5ErrorHandler &) = delete;
  ScopedH5ErrorHandler &operator=(const ScopedH5ErrorHandler &) = delete;

 private:
  H5E_auto2_t func_ = nullptr;
  void *data_ = nullptr;
};

inline void FailH5(const std::string &what) {
  const std::string message =
      "photons multigroup: HDF5 error while " + what + "\n";
  OPAC_ERROR(message.c_str());
}

inline void RequireH5Success(const herr_t status, const std::string &what) {
  if (status != H5_SUCCESS) FailH5(what);
}

inline hid_t OpenFileRead(const std::string &filename) {
  if (!std::filesystem::exists(filename)) {
    const std::string message =
        "photons multigroup: HDF5 file does not exist: " + filename + "\n";
    OPAC_ERROR(message.c_str());
  }
  const hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0) FailH5("opening " + filename + " for reading");
  return file;
}

inline hid_t OpenFileWrite(const std::string &filename, const bool append) {
  const bool reuse = append && std::filesystem::exists(filename);
  const hid_t file = reuse
                         ? H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT)
                         : H5Fcreate(filename.c_str(), H5F_ACC_TRUNC,
                                     H5P_DEFAULT, H5P_DEFAULT);
  if (file < 0) {
    FailH5(std::string(reuse ? "opening " : "creating ") + filename +
           " for writing");
  }
  return file;
}

inline void CloseFile(const hid_t file) {
  RequireH5Success(H5Fclose(file), "closing HDF5 file");
}

inline void CloseGroup(const hid_t group) {
  RequireH5Success(H5Gclose(group), "closing HDF5 group");
}

inline hid_t CreateOrOpenGroup(const hid_t file, const std::string &path,
                               const bool reuse_existing) {
  const bool reuse =
      reuse_existing && H5Lexists(file, path.c_str(), H5P_DEFAULT) > 0;
  const hid_t group = reuse ? H5Gopen(file, path.c_str(), H5P_DEFAULT)
                            : H5Gcreate(file, path.c_str(), H5P_DEFAULT,
                                        H5P_DEFAULT, H5P_DEFAULT);
  if (group < 0) {
    FailH5(std::string(reuse ? "opening" : "creating") + " material group " +
           path);
  }
  return group;
}

inline void SetAttribute(const hid_t file, const std::string &path,
                         const char *name, const int value) {
  RequireH5Success(H5LTset_attribute_int(file, path.c_str(), name, &value, 1),
                   "writing attribute " + std::string(name) + " on " + path);
}

inline void SetAttribute(const hid_t file, const std::string &path,
                         const char *name, const std::string &value) {
  RequireH5Success(
      H5LTset_attribute_string(file, path.c_str(), name, value.c_str()),
      "writing attribute " + std::string(name) + " on " + path);
}

// State for the root-group scan that matches an opacid attribute. A group whose
// own name is the opacid is tracked separately, since it claims that opacid
// whether or not it also carries the attribute.
struct MaterialOpacidSearch {
  int requested = 0;
  std::string direct_name;
  std::string path;
  int matches = 0;
  bool direct_matched = false;
};

inline bool LinkIsGroup(const hid_t loc, const char *link_name) {
  hid_t candidate = H5Gopen(loc, link_name, H5P_DEFAULT);
  if (candidate < 0) return false;
  H5Gclose(candidate);
  return true;
}

inline herr_t FindMaterialByOpacid(hid_t file, const char *link_name,
                                   const H5L_info_t *, void *opaque) {
  auto *search = static_cast<MaterialOpacidSearch *>(opaque);
  if (!LinkIsGroup(file, link_name) ||
      H5Aexists_by_name(file, link_name, SP5::Material::opacid, H5P_DEFAULT) <=
          0) {
    return 0;
  }

  int candidate = 0;
  if (H5LTget_attribute_int(file, link_name, SP5::Material::opacid,
                            &candidate) < 0 ||
      candidate != search->requested) {
    return 0;
  }
  search->matches += 1;
  if (link_name == search->direct_name) search->direct_matched = true;
  if (search->path.empty()) search->path = "/" + std::string(link_name);
  // Keep iterating so that duplicate opacids are detected rather than masked.
  return 0;
}

// Opens a group that H5Lexists has already reported as present, so a failure
// here is a real HDF5 error rather than a missing material.
inline hid_t OpenExistingGroup(const hid_t file, const std::string &path) {
  const hid_t group = H5Gopen(file, path.c_str(), H5P_DEFAULT);
  if (group < 0) FailH5("opening material group " + path);
  return group;
}

inline hid_t OpenMaterialGroupByName(const hid_t file,
                                     const std::string &name) {
  ScopedH5ErrorHandler h5_errors;
  const std::string path = "/" + name;
  if (H5Lexists(file, path.c_str(), H5P_DEFAULT) <= 0 ||
      !LinkIsGroup(file, path.c_str())) {
    return MaterialNotFound;
  }
  return OpenExistingGroup(file, path);
}

inline hid_t OpenMaterialGroupByOpacid(const hid_t file, const int opacid) {
  ScopedH5ErrorHandler h5_errors;
  const std::string direct_name = std::to_string(opacid);
  const std::string direct_path = "/" + direct_name;
  const bool direct = H5Lexists(file, direct_path.c_str(), H5P_DEFAULT) > 0 &&
                      LinkIsGroup(file, direct_path.c_str());

  // Scanning happens even when /<opacid> exists, so that a group named for the
  // opacid and a differently named group carrying it as an attribute collide.
  MaterialOpacidSearch search;
  search.requested = opacid;
  search.direct_name = direct_name;
  const hid_t root = H5Gopen(file, "/", H5P_DEFAULT);
  if (root < 0) FailH5("opening the root group");
  hsize_t index = 0;
  const herr_t status = H5Literate(root, H5_INDEX_NAME, H5_ITER_NATIVE, &index,
                                   FindMaterialByOpacid, &search);
  RequireH5Success(H5Gclose(root), "closing the root group");
  if (status < 0) FailH5("scanning materials for opacid " + direct_name);

  const int claimants =
      search.matches + (direct && !search.direct_matched ? 1 : 0);
  if (claimants > 1) return MaterialAmbiguous;
  if (direct) return OpenExistingGroup(file, direct_path);
  return search.path.empty() ? MaterialNotFound
                             : OpenExistingGroup(file, search.path);
}

struct MaterialSelector {
  std::optional<int> opacid;
  std::string name;
};

inline hid_t OpenMaterialGroup(const hid_t file,
                               const MaterialSelector &selector) {
  if (!selector.name.empty()) {
    return OpenMaterialGroupByName(file, selector.name);
  }
  if (selector.opacid.has_value()) {
    return OpenMaterialGroupByOpacid(file, *selector.opacid);
  }
  return MaterialNotFound;
}

inline void SaveDataBox(const hid_t material, const char *field,
                        const MeanUtilsDataBox &data) {
  RequireH5Success(data.saveHDF(material, field),
                   "writing " + std::string(field));
}

inline void LoadDataBox(const hid_t material, const char *field,
                        MeanUtilsDataBox &data) {
  RequireH5Success(data.loadHDF(material, field),
                   "reading " + std::string(field));
}
#endif

// Log/anti-log transforms used to store and interpolate opacities. A small
// floor keeps toLog well-defined at zero.
PORTABLE_INLINE_FUNCTION
Real ToLog(const Real x) {
  constexpr Real eps = 10.0 * std::numeric_limits<Real>::min();
  return std::log10(std::abs(x) + eps);
}

PORTABLE_INLINE_FUNCTION
Real FromLog(const Real lx) { return std::pow(10., lx); }

// Group-bound access, supporting both a generic indexer (used at construction
// time) and a stored DataBox (used after the bounds are cached).
template <typename GroupBoundsIndexer>
PORTABLE_INLINE_FUNCTION Real
GroupBoundAt(const GroupBoundsIndexer &group_bounds, const int group) {
  return group_bounds[group];
}

PORTABLE_INLINE_FUNCTION
Real GroupBoundAt(const MeanUtilsDataBox &group_bounds, const int group) {
  return group_bounds(group);
}

// Locate the group index containing nu via binary search over the cached
// bounds. Boundary cases are shortcut: nu at or below the first bound maps to
// group 0, and nu at or above the final bound maps to the last group. Callers
// are responsible for range-checking nu before invoking this helper.
PORTABLE_INLINE_FUNCTION
int GroupOfNuImpl(const MeanUtilsDataBox &groupBounds, const int ngroups,
                  const Real nu) {
  // Shortcuts for boundary cases
  if (nu <= GroupBoundAt(groupBounds, 0)) {
    return 0;
  }
  if (nu >= GroupBoundAt(groupBounds, ngroups)) {
    return ngroups - 1;
  }
  // Binary search to find group index containing nu
  int lower = 0;
  int upper = ngroups;
  while (upper - lower > 1) {
    const int middle = (lower + upper) / 2;
    if (nu < GroupBoundAt(groupBounds, middle)) {
      upper = middle;
    } else {
      lower = middle;
    }
  }
  return lower;
}

// Copy user-provided group bounds into the class-owned DataBox.
template <typename GroupBoundsIndexer>
void SetGroupBounds(MeanUtilsDataBox &groupBounds,
                    const GroupBoundsIndexer &group_bounds, const int ngroups) {
  groupBounds.resize(ngroups + 1);
  for (int group = 0; group <= ngroups; ++group) {
    groupBounds(group) = GroupBoundAt(group_bounds, group);
  }
}

// Copy the class-owned group bounds back out (e.g. for HDF5 export).
inline void ExportGroupBounds(MeanUtilsDataBox &groupBounds,
                              const MeanUtilsDataBox &storedBounds,
                              const int ngroups) {
  groupBounds.resize(ngroups + 1);
  for (int group = 0; group <= ngroups; ++group) {
    groupBounds(group) = storedBounds(group);
  }
}

// Validate half-open group bounds [nu_g, nu_{g+1}): strictly increasing,
// nonnegative, with only the final bound permitted to be IEEE +infinity.
template <typename GroupBoundsIndexer>
void ValidateGroupBounds(const GroupBoundsIndexer &group_bounds,
                         const int ngroups) {
  if (ngroups <= 0) {
    OPAC_ERROR("photons multigroup: ngroups must be positive");
  }
  for (int group = 0; group <= ngroups; ++group) {
    const Real bound = GroupBoundAt(group_bounds, group);
    if (std::isnan(bound)) {
      OPAC_ERROR("photons multigroup: group bounds must be finite "
                 "or IEEE +infinity");
    }
    if (std::isinf(bound) && bound < 0.) {
      OPAC_ERROR("photons multigroup: group bounds may not be -infinity");
    }
    if (group == 0) {
      if (!(bound >= 0.)) {
        OPAC_ERROR("photons multigroup: first group bound must be "
                   "nonnegative");
      }
    } else if (!(bound > GroupBoundAt(group_bounds, group - 1))) {
      OPAC_ERROR("photons multigroup: group bounds must be strictly "
                 "increasing");
    }
    if (!std::isfinite(bound) && group != ngroups) {
      OPAC_ERROR("photons multigroup: only the final group bound may be "
                 "IEEE +infinity");
    }
  }
}

// Sample a group's frequency range on a logarithmic, midpoint-rule grid,
// invoking sample_op(nu, weight) for each of NNuPerGroup points. Extreme
// bounds (nu=0, +infinity, or far from the thermal peak) are clamped to a
// thermal-aware window so the integral stays well-conditioned.
template <typename PC, typename SampleOp>
void ForEachGroupFrequencySample(const Real temp, const Real nuMin,
                                 const Real nuMax, const int NNuPerGroup,
                                 SampleOp &&sample_op) {
  // For [0, ∞) or very wide ranges, use thermal-aware sampling
  // For reasonable finite ranges, integrate over the full group bounds
  const Real nu_thermal_min = 1.e-3 * PC::kb * temp / PC::h;
  const Real nu_thermal_max = 1.e3 * PC::kb * temp / PC::h;

  // Determine if we need special handling
  const bool is_lower_extreme = (nuMin == 0.) || (nuMin < 0.1 * nu_thermal_min);
  const bool is_upper_extreme =
      !std::isfinite(nuMax) || (nuMax > 10. * nu_thermal_max);

  // Set integration bounds, but ensure they're valid
  Real nu_sample_min = is_lower_extreme ? nu_thermal_min : nuMin;
  Real nu_sample_max = is_upper_extreme ? nu_thermal_max : nuMax;

  // If thermal-aware bounds are invalid, use a small but valid range within
  // group bounds
  if (nu_sample_min >= nu_sample_max) {
    if (std::isfinite(nuMax) && nuMax > 0.) {
      // Group is [0 or small, nuMax]: sample near nuMax
      nu_sample_min = 0.5 * nuMax;
      nu_sample_max = nuMax;
    } else {
      // Group extends to infinity: sample around thermal peak
      nu_sample_min = 0.1 * nu_thermal_max;
      nu_sample_max = nu_thermal_max;
    }
  }

  // Use logarithmic spacing with midpoint rule
  const Real lNuMin = ToLog(nu_sample_min);
  const Real lNuMax = ToLog(nu_sample_max);
  const Real dlnu = (lNuMax - lNuMin) / NNuPerGroup;
  for (int inu = 0; inu < NNuPerGroup; ++inu) {
    const Real lnu = lNuMin + (inu + 0.5) * dlnu;
    const Real nu = FromLog(lnu);
    sample_op(nu, nu * dlnu);
  }
}

// Evaluate the Planck function B_nu and its temperature derivative dB_nu/dT,
// switching to the Wien closed form deep in the tail (see wien_tail_x).
template <typename PC>
void ThermalWeightsAtNu(const PlanckDistribution<PC> &dist, const Real temp,
                        const Real nu, Real &B, Real &dBdT) {
  const Real x = PC::h * nu / (PC::kb * temp);
  if (x < wien_tail_x) {
    B = dist.ThermalDistributionOfTNu(temp, nu);
    dBdT = dist.DThermalDistributionOfTNuDT(temp, nu);
    return;
  }

  const Real expMinusX = std::exp(-x);
  B = (2. * PC::h * nu * nu * nu / (PC::c * PC::c)) * expMinusX;
  dBdT = 2. * PC::h * PC::h * nu * nu * nu * nu * expMinusX /
         (temp * temp * PC::c * PC::c * PC::kb);
}

} // namespace impl
} // namespace photons
} // namespace singularity

#endif // SINGULARITY_OPAC_PHOTONS_MEAN_PHOTON_UTILS_
