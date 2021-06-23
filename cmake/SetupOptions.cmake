# Options

include(CMakeDependentOption)

#=======================================
# Build options
#=======================================
option (SINGULARITY_INSTALL_LIBRARY "Enable installing the libraries to default locations" OFF)
option (SINGULARITY_USE_FORTRAN "Enable fortran bindings" OFF)
option (SINGULARITY_HIDE_MORE_WARNINGS "hide more warnings" OFF)
option (SINGULARITY_BUILD_TESTS "Compile tests" OFF)
option (SINGULARITY_BETTER_DEBUG_FLAGS
  "Better debug flags for singularity" ON)

#=======================================
# Dependency options
#=======================================
# check for using in-tree third-party dependencies
option (SINULARITY_KOKKOS_IN_TREE "Use in-tree dependencies" OFF)

option (SINGULARITY_USE_KOKKOS "Use Kokkos for portability" OFF)
option (SINGULARITY_USE_CUDA "Enable cuda support" OFF)
option (SINGULARITY_USE_HDF5 "Pull in hdf5" OFF)

cmake_dependent_option(SINGULARITY_USE_MPI "Link to MPI" OFF "HDF5_IS_PARALLEL" ON)

#=======================================
# Options logic
#=======================================
# check for currently incompatible option
if(SINGULARITY_USE_CUDA)
  if(NOT CUDAToolkit_FOUND)
    message(FATAL_ERROR "Requested build with CUDA, but could not find in the environment.")
  endif()
  if(SINGULARITY_USE_KOKKOS)
    if(NOT Kokkos_FOUND)
      message(WARNING "Requested build with Kokkos, but could not find in the environment.")
      message(STATUS  "Configuring internal Kokkos submodule.")
      set(NVCC_WRAPPER_BIN ${CMAKE_CURRENT_SOURCE_DIR}/utils/kokkos/bin/nvcc_wrapper CACHE STRING "")
      set(Kokkos_ENABLE_CUDA ON CACHE BOOL "" FORCE)
      set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "" FORCE)
      set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "" FORCE)
      set(Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE ON CACHE BOOL "" FORCE)
      add_subdirectory(${PROJECT_SOURCE_DIR}/utils/kokkos)
    endif()
  else()
    message(FATAL_ERROR "Cuda without kokkos is not currently supported")
  endif()
endif()
