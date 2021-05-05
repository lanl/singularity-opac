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

add_library(catch2_define catch2_define.cpp)
target_link_libraries(catch2_define PUBLIC
                      Catch2::Catch2
                      singularity-opac::libs singularity-opac::flags)

list(APPEND unit_tests_SOURCES
  test_trivial.cpp
  test_gray_opacities.cpp
)

add_executable(unit_tests "${unit_tests_SOURCES}")
target_link_libraries(unit_tests PRIVATE catch2_define
                      Catch2::Catch2
                      singularity-opac::libs singularity-opac::flags)

include(Catch)
catch_discover_tests(unit_tests)
