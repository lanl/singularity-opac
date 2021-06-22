#=======================================
# Compiler logic
#=======================================

#TODO: Do CMake CUDA options help here?
#if(SINGULARITY_USE_KOKKOS AND SINGULARITY_USE_CUDA)
# set nvcc wrapper default compiler
#  if(NOT "$ENV{CXX}x" STREQUAL "x" AND
#     "$ENV{NVCC_WRAPPER_DEFAULT_COMPILER}x" STREQUAL "x")
#    set(ENV{NVCC_WRAPPER_DEFAULT_COMPILER} "$ENV{CXX}")
#  endif()
  # set necessary kokkos build options if building inline
#  if(NOT Kokkos_FOUND)
#    set(NVCC_WRAPPER_DIR ${CMAKE_CURRENT_SOURCE_DIR}/utils/kokkos/bin CACHE STRING "")
#  endif()

 # set(CMAKE_CXX_COMPILER ${KOKKOS_CXX_COMPILER} CACHE STRING "")
#endif()

# set C++ options
enable_language(CXX)
include(CMakeDetermineCXXCompiler)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Fortran options
if(SINGULARITY_USE_FORTRAN)
  enable_language(Fortran)
  include(CMakeDetermineFortranCompiler) #NOTE: this can be annoying with Cray machines
endif()

