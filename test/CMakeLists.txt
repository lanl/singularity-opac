# © 2021. Triad National Security, LLC. All rights reserved.  This
# program was produced under U.S. Government contract 89233218CNA000001
# for Los Alamos National Laboratory (LANL), which is operated by Triad
# National Security, LLC for the U.S.  Department of Energy/National
# Nuclear Security Administration. All rights in the program are
# reserved by Triad National Security, LLC, and the U.S. Department of
# Energy/National Nuclear Security Administration. The Government is
# granted for itself and others acting on its behalf a nonexclusive,
# paid-up, irrevocable worldwide license in this material to reproduce,
# prepare derivative works, distribute copies to the public, perform
# publicly and display publicly, and to permit others to do so.

if(NOT Catch2_FOUND)
  # Get catch2
  message(STATUS "Fetching Catch2 as needed")
  Include(FetchContent)
  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v2.13.10)
  FetchContent_MakeAvailable(Catch2)
  list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/contrib)
endif()

# Build tests
message(STATUS "Configuring unit tests")
#find_package(${PROJECT_NAME} CONFIG REQUIRED)

add_library(catch2_define)
target_sources(catch2_define PUBLIC catch2_define.cpp)
target_link_libraries(catch2_define PUBLIC
                      Catch2::Catch2)

add_executable(${PROJECT_NAME}_unit_tests)
target_sources(${PROJECT_NAME}_unit_tests
PRIVATE
  test_thermal_dist.cpp
  test_scalefree_opacities.cpp
  test_gray_opacities.cpp
  test_epbremsstrahlung_opacities.cpp
  test_gray_s_opacities.cpp
  test_brt_opacities.cpp
  test_thomson_s_opacities.cpp
  test_chebyshev.cpp
  test_spiner_opac_neutrinos.cpp
  test_mean_opacities.cpp
  test_variant.cpp
)

target_link_libraries(${PROJECT_NAME}_unit_tests
PRIVATE
  catch2_define
  ${PROJECT_NAME})

# Ensure code works with C++11 and earlier
# TODO(MM): Remove this later when it's not needed.
set_target_properties(${PROJECT_NAME}_unit_tests
  PROPERTIES CXX_STANDARD 14
  CXX_STANDARD_REQUIRED YES
  CXX_EXTENSIONS NO)

include(Catch)
catch_discover_tests(${PROJECT_NAME}_unit_tests)
