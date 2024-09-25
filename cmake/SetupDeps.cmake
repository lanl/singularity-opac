include(FeatureSummary)
include(CMakeDependentOption)
#=======================================
# Setup CUDAToolkit
# - provideds CUDA::toolkit
#=======================================
if (NOT TARGET CUDA::toolkit)
  find_package(CUDAToolkit QUIET)
else()
  message(status "CUDA::toolkit provided by parent package")
endif()

#=======================================
# Setup ports of call
# - provides PortsofCall::PortsofCall
#=======================================
find_package(PortsofCall REQUIRED)
target_link_libraries(singularity-opac::flags INTERFACE PortsofCall::PortsofCall)

#=======================================
# Setup Kokkos
# - provides Kokkos::kokkos
#=======================================
if (SINGULARITY_USE_KOKKOS)
  if (NOT TARGET Kokkos::kokkos)
    add_subdirectory(utils/kokkos)
    find_package(Kokkos QUIET)
  else()
    message(status "Kokkos::kokkos provided by parent package")
  endif()
endif()

#=======================================
# Find HDF5
# - cmake@3.20+ provides HDF5::HDF5, but
#   prior versions do not
#=======================================
# Can't play the target games above unless using
# cmake 3.20+, so leave this one for now.
# JMM: DO NOT SEARCH FOR HDF5 IF WE DON'T NEED IT
if (SINGULARITY_USE_HDF5)
find_package(HDF5 COMPONENTS C HL)

# findpackage doesnt export an interface for HDF5,
# so create one
if (HDF5_FOUND)
    target_compile_definitions(singularity-opac::flags INTERFACE SPINER_USE_HDF)
    add_library(${PROJECT_NAME}::hdf5 INTERFACE IMPORTED)
    set_target_properties(${PROJECT_NAME}::hdf5
    PROPERTIES
        INTERFACE_LINK_LIBRARIES "${HDF5_LIBRARIES};${HDF5_HL_LIBRARIES}"
        INTERFACE_COMPILE_DEFINITIONS "SINGULARITY_USE_HDF5;SPINER_USE_HDF"
        INTERFACE_INCLUDE_DIRECTORIES "${HDF5_INCLUDE_DIRS}"
    )

    # if HDF5 is parallel, also get MPI libraries
    if (HDF5_IS_PARALLEL)
      message(status "Parallel HDF5 found. Looking for MPI")
      # Provides MPI::MPI_CXX
      if (NOT TARGET MPI::MPI_CXX)
	find_package(MPI COMPONENTS)# CXX REQUIRED)
	if (NOT MPI_FOUND)
	  message(FATAL_ERROR
	    "HDF5 parallel found, but could not find MPI! "
	    "Try specifying a path to an MPI directory or MPI compiler "
	    " via -DMPI_HOME=/path or -DMPI_CXX_COMPILER=/path, "
	    "or set CMAKE_PREFIX_PATH to your preferred HDF5 directory "
	    "via -DCMAKE_PREFIX_PATH=/path.")
	endif()
      else()
	message("MPI::MPI_CXX provided by parent package")
      endif()
    else()
      message(status "HDF5 is not parallel")
    endif()
else()
  message(STATUS "No HDF5")
endif()
endif()

# JMM: I'm putting this here, as it depends on what happens with HDF5,
# and I want that to be set AFTER we set all our options, so we don't
# set HDF5 unnecessarily.
cmake_dependent_option(SINGULARITY_USE_MPI
        "Link to MPI"
        ON "${SINGULARITY_USE_MPI};${HDF5_IS_PARALLEL}"
        OFF)

#=======================================
# Setup Catch2
# - provides Catch2::Catch2
#=======================================
if (NOT TARGET Catch2::Catch2)
  find_package(Catch2 QUIET)
endif()
