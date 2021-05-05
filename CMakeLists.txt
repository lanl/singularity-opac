# © 2021. Triad National Security, LLC. All rights reserved.  This
# program was produced under U.S. Government contract
# 89233218CNA000001 for Los Alamos National Laboratory (LANL), which
# is operated by Triad National Security, LLC for the U.S.  Department
# of Energy/National Nuclear Security Administration. All rights in
# the program are reserved by Triad National Security, LLC, and the
# U.S. Department of Energy/National Nuclear Security
# Administration. The Government is granted for itself and others
# acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
# license in this material to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.

cmake_minimum_required(VERSION 3.14)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
set(CMAKE_EXPORT_COMPILE_COMMANDS On)

# Options
option (SINGULARITY_USE_KOKKOS "Use Kokkos for portability" OFF)
option (SINGULARITY_USE_FORTRAN "Enable fortran bindings" OFF)
option (SINGULARITY_USE_CUDA "Enable cuda support" OFF)
option (SINGULARITY_USE_HDF5 "Pull in hdf5" OFF)
option (SINGULARITY_HIDE_MORE_WARNINGS "hide more warnings" OFF)
option (SINGULARITY_BUILD_TESTS "Compile tests" OFF)
option (SINGULARITY_BETTER_DEBUG_FLAGS
  "Better debug flags for singularity" ON)

# Don't allow in-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
  message(FATAL_ERROR
    "You cannot build in a source directory (or any directory with a CMakeLists.txt file). "
    "Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

# if cuda, need to set cxx compiler to nvcc_wrapper before compilers
# are checked
if(SINGULARITY_USE_KOKKOS AND SINGULARITY_USE_CUDA)
# set nvcc wrapper default compiler
  if(NOT "$ENV{CXX}x" STREQUAL "x" AND
     "$ENV{NVCC_WRAPPER_DEFAULT_COMPILER}x" STREQUAL "x")
    set(ENV{NVCC_WRAPPER_DEFAULT_COMPILER} "$ENV{CXX}")
  endif()
  # set necessary kokkos build options if building inline
  if(NOT SINGULARITY_KOKKOS_INSTALL_DIR)
    set(NVCC_WRAPPER_DIR ${CMAKE_CURRENT_SOURCE_DIR}/utils/kokkos/bin CACHE STRING "")
  else()
    set(NVCC_WRAPPER_DIR ${SINGULARITY_KOKKOS_INSTALL_DIR}/bin CACHE STRING "")
  endif()
  set(CMAKE_CXX_COMPILER ${NVCC_WRAPPER_DIR}/nvcc_wrapper CACHE STRING "")
endif()

project(singularity-opac VERSION 0.0.1)

if(SINGULARITY_USE_FORTRAN)
  enable_language(CXX Fortran)
else()
  enable_language(CXX)
endif()
include(CMakeDetermineCXXCompiler)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Patches variant to be compatible with cuda
# Assumes "patch" is present on system
message(STATUS "Patching mpark::variant to support GPUs")
execute_process(COMMAND patch -N -s -V never
                ${CMAKE_CURRENT_SOURCE_DIR}/utils/variant/include/mpark/variant.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/utils/cuda_compatibility.patch)

add_library (singularity-opac::flags INTERFACE IMPORTED GLOBAL)
add_library (singularity-opac::libs INTERFACE IMPORTED)

# If the user doesn't specify a build type, prefer RelWithDebInfo
set(default_build_type "RelWithDebInfo")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
    STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# xl fix
target_compile_options(singularity-opac::flags INTERFACE
                       $<$<COMPILE_LANG_AND_ID:CXX,XL>:-std=c++1y;-qxflag=disable__cplusplusOverride>)
target_link_options(singularity-opac::flags INTERFACE
                   $<$<COMPILE_LANG_AND_ID:CXX,XL>:-std=c++1y;-qxflag=disable__cplusplusOverride>)

# Base Include directories
target_include_directories(singularity-opac::libs
  INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/utils>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/opacities>)
target_include_directories(singularity-opac::flags
  INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/utils>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/opacities>)

# check for currently incompatible option
if(SINGULARITY_USE_CUDA AND NOT SINGULARITY_USE_KOKKOS)
  message(FATAL_ERROR "Cuda without kokkos is not currently supported")
endif()

if (SINGULARITY_USE_HDF5)
  target_compile_definitions(singularity-opac::flags INTERFACE
                             SPINER_USE_HDF)
  if(SINGULARITY_HDF5_INSTALL_DIR)
    list(APPEND CMAKE_PREFIX_PATH "${SINGULARITY_HDF5_INSTALL_DIR}")
  endif()
  find_package(HDF5 COMPONENTS C HL)
  if(HDF5_FOUND)
    add_library(singularity-opac::hdf5 INTERFACE IMPORTED GLOBAL)
    set_target_properties(singularity-opac::hdf5 PROPERTIES
      INTERFACE_LINK_LIBRARIES "${HDF5_LIBRARIES};${HDF5_HL_LIBRARIES}"
      INTERFACE_COMPILE_DEFINITIONS "SINGULARITY_USE_HDF5"
      INTERFACE_INCLUDE_DIRECTORIES "${HDF5_INCLUDE_DIRS}")
    if(HDF5_IS_PARALLEL)
      if(SINGULARITY_MPI_INSTALL_DIR)
        list(APPEND CMAKE_PREFIX_PATH "${SINGULARITY_MPI_INSTALL_DIR}")
      endif()
      find_package(MPI COMPONENTS CXX)
      if(MPI_FOUND)
        target_include_directories(singularity-opac::libs INTERFACE "${MPI_CXX_INCLUDE_DIRS}")
      endif()
    endif()
  else()
    message(FATAL_ERROR "HDF5 was requested but not found. Can be disabled with -DSINGULARITY_USE_HDF5=OFF")
  endif()
  target_link_libraries(singularity-opac::libs INTERFACE singularity-opac::hdf5)
endif()

if (SINGULARITY_USE_KOKKOS)
  target_compile_definitions(singularity-opac::flags INTERFACE
                             PORTABILITY_STRATEGY_KOKKOS)
  # Import Kokkos if not already available as a build target
  if (NOT TARGET Kokkos::kokkos)
    if (SINGULARITY_KOKKOS_INSTALL_DIR)
      list(APPEND CMAKE_PREFIX_PATH "${SINGULARITY_KOKKOS_INSTALL_DIR}")
      set(Kokkos_ROOT "${SINGULARITY_KOKKOS_INSTALL_DIR}/lib64/cmake/Kokkos")
      find_package(Kokkos REQUIRED)
    else()
      if(SINGULARITY_USE_CUDA)
        set(Kokkos_ENABLE_CUDA ON CACHE BOOL "" FORCE)
        set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "" FORCE)
        set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "" FORCE)
        set(Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE ON CACHE BOOL "" FORCE)
      endif()
      add_subdirectory(${PROJECT_SOURCE_DIR}/utils/kokkos)
      #set(Kokkos_ROOT ${PROJECT_BINARY_DIR}/utils/kokkos)
    endif()
  endif()
  # cuda specific flags
  if(SINGULARITY_USE_CUDA)
    if(SINGULARITY_HIDE_MORE_WARNINGS)
      target_compile_options(
        singularity-opac::flags
         INTERFACE # Generator expression shamelessly copied from EAP
         "$<$<COMPILE_LANGUAGE:CXX>:--expt-relaxed-constexpr;>"
           )
    else()
      target_compile_options(
        singularity-opac::flags
         INTERFACE # Generator expression shamelessly copied from EAP
         "$<$<COMPILE_LANGUAGE:CXX>:--expt-relaxed-constexpr;-Xcudafe;--diag_suppress=esa_on_defaulted_function_ignored;>"
           )
    endif()
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
      target_compile_options(singularity-opac::flags INTERFACE
                             -use_fast_math)
    elseif(CMAKE_BUILD_TYPE STREQUAL "Debug" AND SINGULARITY_BETTER_DEBUG_FLAGS)
      target_compile_options(
        singularity-opac::flags
        INTERFACE
        "$<$<COMPILE_LANGUAGE:CXX>:-G;-lineinfo;>"
        )
    endif()
  endif()
  target_link_libraries(singularity-opac::libs INTERFACE Kokkos::kokkos)
else()
  target_include_directories(singularity-opac::libs INTERFACE
                             $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/utils/herumi-fmath>)
endif ()

# TODO(JMM)
# Add source

# Unit tests
if (SINGULARITY_BUILD_TESTS)

  include(CTest)

  # Get catch2
  message(STATUS "Fetching Catch2 as needed")
  Include(FetchContent)
  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v2.13.1)
  FetchContent_MakeAvailable(Catch2)
  list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/contrib)

  # Build tests
  message(STATUS "Configuring unit tests")
  add_subdirectory(test)
endif()

# TODO(JMM)
# Add install commands
