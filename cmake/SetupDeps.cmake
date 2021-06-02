
#=======================================
# Find Kokkos
#=======================================
find_package(Kokkos QUIET)

#=======================================
# Find CUDA
#=======================================
#find_package(Cuda)

#=======================================
# Find HDF5
#=======================================
find_package(HDF5 COMPONENTS C HL QUIET)

# findpackage doesnt export an interface for HDF5,
# so create one
if (HDF5_FOUND)
    add_library( singularity-opac::hdf5 INTERFACE IMPORTED)
    set_target_properties(singularity-opac::hdf5
    PROPERTIES
        INTERFACE_LINK_LIBRARIES "${HDF5_LIBRARIES};${HDF5_HL_LIBRARIES}"
        INTERFACE_COMPILE_DEFINITIONS "SINGULARITY_USE_HDF5"
        INTERFACE_INCLUDE_DIRECTORIES "${HDF5_INCLUDE_DIRS}"
    )

endif()

#=======================================
# Find MPI
#=======================================
find_package(MPI COMPONENTS CXX QUIET)
