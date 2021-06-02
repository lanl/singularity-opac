# Options

include (CMakeDependentOption)

#=======================================
# Build options
#=======================================
option (SINGULARITY_USE_FORTRAN "Enable fortran bindings" OFF)
option (SINGULARITY_HIDE_MORE_WARNINGS "hide more warnings" OFF)
option (SINGULARITY_BUILD_TESTS "Compile tests" OFF)
option (SINGULARITY_BETTER_DEBUG_FLAGS
  "Better debug flags for singularity" ON)

#=======================================
# Dependency options
#=======================================
cmake_dependent_option (SINGULARITY_USE_KOKKOS "Use Kokkos for portability" OFF "Kokkos_FOUND" ON)
cmake_dependent_option (SINGULARITY_USE_CUDA "Enable cuda support" OFF "CUDA_FOUND" ON)
cmake_dependent_option (SINGULARITY_USE_HDF5 "Pull in hdf5" OFF "HDF5_FOUND" ON)

#=======================================
# Options logic
#=======================================
# check for currently incompatible option
if(SINGULARITY_USE_CUDA AND NOT SINGULARITY_USE_KOKKOS)
  message(FATAL_ERROR "Cuda without kokkos is not currently supported")
endif()