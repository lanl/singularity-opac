include(FeatureSummary)
#=======================================
# Setup CUDAToolkit
# - provideds CUDA::toolkit
#=======================================
find_package(CUDAToolkit QUIET)

#=======================================
# Setup Kokkos
# - provides Kokkos::kokkos
#=======================================
find_package(Kokkos QUIET)

#=======================================
# Find HDF5
# - cmake@3.20+ provides HDF5::HDF5, but
#   prior versions do not
#=======================================
find_package(HDF5 COMPONENTS C HL QUIET)

# findpackage doesnt export an interface for HDF5,
# so create one
if (HDF5_FOUND)
    add_library(${PROJECT_NAME}::hdf5 INTERFACE IMPORTED)
    set_target_properties(${PROJECT_NAME}::hdf5
    PROPERTIES
        INTERFACE_LINK_LIBRARIES "${HDF5_LIBRARIES};${HDF5_HL_LIBRARIES}"
        INTERFACE_COMPILE_DEFINITIONS "SINGULARITY_USE_HDF5;SPINER_USE_HDF"
        INTERFACE_INCLUDE_DIRECTORIES "${HDF5_INCLUDE_DIRS}"
    )

    # if HDF5 is parallel, also get MPI libraries
    if (HDF5_IS_PARALELL)
        find_package(MPI COMPONENTS CXX REQUIRED)
    endif()
endif()

#=======================================
# Setup Catch2
# - provides Catch2::Catch2
#=======================================
find_package(Catch2 QUIET)



