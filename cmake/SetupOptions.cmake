# Options

include (CMakeDependentOption)

#=======================================
# Build options
#=======================================
option (SINGULARITY_INSTALL_LIBRARY "Enable installing the libraries to default locations")
option (SINGULARITY_USE_FORTRAN "Enable fortran bindings" OFF)
option (SINGULARITY_HIDE_MORE_WARNINGS "hide more warnings" OFF)
option (SINGULARITY_BUILD_TESTS "Compile tests" OFF)
option (SINGULARITY_BETTER_DEBUG_FLAGS
  "Better debug flags for singularity" ON)

#=======================================
# Dependency options
#=======================================
# check for using in-tree third-party dependencies
option (SINULARITY_IN_TREE "Use in-tree dependencies" ON)

cmake_dependent_option (SINGULARITY_USE_KOKKOS "Use Kokkos for portability" ON "Kokkos_FOUND" OFF)
cmake_dependent_option (SINGULARITY_USE_CUDA "Enable cuda support" ON "SINGULARITY_USE_KOKKOS" OFF)
cmake_dependent_option (SINGULARITY_USE_HDF5 "Pull in hdf5" ON "HDF5_FOUND" OFF)

#=======================================
# Options logic
#=======================================
# check for currently incompatible option
if(SINGULARITY_USE_CUDA AND NOT SINGULARITY_USE_KOKKOS)
  message(FATAL_ERROR "Cuda without kokkos is not currently supported")
endif()