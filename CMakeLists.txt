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

project("singularity-opac"
      VERSION 0.0.1
      DESCRIPTION
        ""
      HOMEPAGE_URL "https://github.com/lanl/singularity-opac/"
)

# Don't allow in-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
  message(FATAL_ERROR
    "You cannot build in a source directory (or any directory with a CMakeLists.txt file). "
    "Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

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

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

add_library(${PROJECT_NAME} INTERFACE IMPORTED GLOBAL)
# so we can use with add_subdirectory
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

include (SetupDeps)
include (SetupOptions)
include (SetupCompilers)
include (SetupFlags)

include(GNUInstallDirs)
include(CTest)

# flag if this is the main project or being included
set(IS_TOPLEVEL_PROJECT FALSE)
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  set(IS_TOPLEVEL_PROJECT TRUE)
endif()


set(CMAKE_EXPORT_COMPILE_COMMANDS On)


# Patches variant to be compatible with cuda
# Assumes "patch" is present on system
message(STATUS "Patching mpark::variant to support GPUs")
execute_process(COMMAND patch -N -s -V never
                ${CMAKE_CURRENT_SOURCE_DIR}/utils/variant/include/mpark/variant.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/utils/cuda_compatibility.patch)

# Flags and libraries for building, see cmake/SetupFlags.cmake
#add_library (singularity-opac::flags INTERFACE IMPORTED)


# Invoke CMake scripts for setup
if(SINGULARITY_INSTALL_LIBRARY)
  include (SetupInstall)
endif()

if(SINGULARITY_BUILD_TESTS)
  add_subdirectory(test)
endif()

set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")

include(CPack)

