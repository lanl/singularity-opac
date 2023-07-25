include(FeatureSummary)
# =======================================
# Setup CUDAToolkit - provideds CUDA::toolkit
# =======================================
if(NOT TARGET CUDA::toolkit)
  find_package(CUDAToolkit QUIET)
else()
  message(status "CUDA::toolkit provided by parent package")
endif()

if(NOT TARGET mpark_variant)
  add_subdirectory(utils/variant)
endif()
target_link_libraries(singularity-opac INTERFACE mpark_variant)
# =======================================
# Setup ports of call - provides PortsofCall::PortsofCall
# =======================================
find_package(ports-of-call QUIET)
if(NOT ports-of-call_FOUND)
  add_subdirectory(utils/ports-of-call)
endif()
target_link_libraries(singularity-opac INTERFACE ports-of-call::ports-of-call)

find_package(spiner QUIET)
if(NOT spiner_FOUND)
  add_subdirectory(utils/spiner)
endif()
target_link_libraries(singularity-opac INTERFACE spiner::spiner)

# =======================================
# Setup Kokkos - provides Kokkos::kokkos
# =======================================
if(NOT TARGET Kokkos::kokkos)
  find_package(Kokkos QUIET)
else()
  message(STATUS "Kokkos::kokkos provided by parent package")
endif()

if(NOT TARGET Kokkos::kokkos OR Kokkos_NOTFOUND)
  message(STATUS "Kokkos not found, will not include in build requirements")
else()
  target_link_libraries(singularity-opac INTERFACE Kokkos::kokkos)
endif()

# =======================================
# Find HDF5 - cmake@3.20+ provides HDF5::HDF5, but prior versions do not
# =======================================
# Can't play the target games above unless using cmake 3.20+, so leave this one
# for now.
find_package(
  HDF5
  COMPONENTS C HL
  QUIET)

# findpackage doesnt export an interface for HDF5, so create one
if(HDF5_FOUND)
  target_compile_definitions(singularity-opac
                             INTERFACE "SINGULARITY_USE_HDF5;SPINER_USE_HDF")
  target_include_directories(singularity-opac INTERFACE "${HDF5_INCLUDE_DIRS}")
  target_link_libraries(singularity-opac
                        INTERFACE "${HDF5_LIBRARIES};${HDF5_HL_LIBRARIES}")

  # if HDF5 is parallel, also get MPI libraries
  if(HDF5_IS_PARALLEL)
    message(status "Parallel HDF5 found. Looking for MPI")
    # Provides MPI::MPI_CXX
    if(NOT TARGET MPI::MPI_CXX)
      find_package(MPI COMPONENTS) # CXX REQUIRED)
      if(NOT MPI_FOUND)
        message(
          FATAL_ERROR
            "HDF5 parallel found, but could not find MPI! "
            "Try specifying a path to an MPI directory or MPI compiler "
            " via -DMPI_HOME=/path or -DMPI_CXX_COMPILER=/path, "
            "or set CMAKE_PREFIX_PATH to your preferred HDF5 directory "
            "via -DCMAKE_PREFIX_PATH=/path.")
      endif()
    else()
      message("MPI::MPI_CXX provided by parent package")
    endif()
    target_link_libraries(singularity-opac INTERFACE MPI::MPI_CXX)
  endif()
endif()

# =======================================
# Setup Catch2 - provides Catch2::Catch2
# =======================================
if(NOT TARGET Catch2::Catch2)
  find_package(Catch2 QUIET)
endif()
