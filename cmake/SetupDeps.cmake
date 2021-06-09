
#=======================================
# Find Kokkos
#=======================================
find_package(Kokkos QUIET)

#=======================================
# Find HDF5
#=======================================
find_package(HDF5 COMPONENTS C HL QUIET)

# findpackage doesnt export an interface for HDF5,
# so create one
if (HDF5_FOUND)
    add_library(${PROJECT_NAME}::hdf5 INTERFACE IMPORTED)
    set_target_properties(${PROJECT_NAME}::hdf5
    PROPERTIES
        INTERFACE_LINK_LIBRARIES "${HDF5_LIBRARIES};${HDF5_HL_LIBRARIES}"
        INTERFACE_COMPILE_DEFINITIONS "SINGULARITY_USE_HDF5"
        INTERFACE_INCLUDE_DIRECTORIES "${HDF5_INCLUDE_DIRS}"
    )

endif()

#=======================================
# Find Catch2
#=======================================
find_package(Catch2 QUIET)

#=======================================
# Find MPI
#=======================================
find_package(MPI COMPONENTS CXX QUIET)
