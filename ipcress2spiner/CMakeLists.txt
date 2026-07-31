#------------------------------------------------------------------------------
# © 2021-2023. Triad National Security, LLC. All rights reserved.  This
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
#------------------------------------------------------------------------------
project(ipcress2spiner VERSION 0.8.0.0 LANGUAGES CXX Fortran)

add_executable(ipcress2spiner
  parse_cli.cpp
  parse_cli.hpp
  make_spiner_databox.hpp
  interpolate_opacity.hpp
  main.cpp
)

# Turn on verbosity right before your package search
#set(CMAKE_FIND_DEBUG_MODE ON)
find_package(gandolf)
#set(CMAKE_FIND_DEBUG_MODE OFF)

target_include_directories(ipcress2spiner
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_BINDIR}>
    ${GANDOLF_INCLUDE_DIRECTORY}
)


target_link_libraries(ipcress2spiner
  PRIVATE
    singularity-opac::singularity-opac
    ${GANDOLF_LIBRARY}
    PortsofCall::PortsofCall
)

set_property(TARGET ipcress2spiner PROPERTY LINKER_LANGUAGE Fortran)
# Define the full version as a string macro
target_compile_definitions(ipcress2spiner PRIVATE
  IPCRESS2SPINER_VERSION=\"${PROJECT_VERSION}\"
)
install(TARGETS ipcress2spiner DESTINATION ${CMAKE_INSTALL_BINDIR})


# TODO: Add tests for sesame2spiner here.
