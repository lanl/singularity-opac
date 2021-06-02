

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

# Easier generator expressions
set(build_debug "$<CONFIG:Debug>")
set(build_release "$<CONFIG:Release>")
set(cxx_lang "$<COMPILE_LANGUAGE:CXX>")
set(cxx_xl "$<COMPILE_LANG_AND_ID:CXX,XL>")
set(with_hdf5 "$<BOOL:SINGULARITY_USE_HDF5>")
set(with_kokkos "$<BOOL:SINGULARITY_USE_KOKKOS>")
set(with_cuda "$<BOOL:SINGULARITY_USE_CUDA>")
set(hide_more_warn "$<BOOL:SINGULARITY_HIDE_MORE_WARNINGS>")
set(better_debug "$<BOOL:SINGULARITY_BETTER_DEBUG_FLAGS>")

# xl fix
target_compile_options(singularity-opac::flags INTERFACE
                        $<${cxx_xl}:
                            "-std=c++1y;-qxflag=disable__cplusplusOverride"
                        >
)
target_link_options(singularity-opac::flags INTERFACE
                        $<${cxx_xl}:
                            "-std=c++1y;-qxflag=disable__cplusplusOverride"
                        >
)

# Base Include directories
target_include_directories(singularity-opac::flags
  INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/utils>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
)

target_compile_options(singularity-opac::flags
INTERFACE
    $<${with_kokkos}:
        $<${with_cuda}:
            $<${cxx_lang}:
                "--expt-relaxed-constexpr"
                $<${hide_more_warn}:
                    "-Xcudafe;--diag_suppress=esa_on_defaulted_function_ignored"
                >
            >
            $<${build_release}:
                "-use_fast_math"
            >
            $<${build_debug}:
                $<${better_debug}:
                    $<${cxx_lang}:
                        "-G;-lineinfo"
                    >
                >
            >
        >
    >
)

# target_link_libraries brings in compile flags, compile defs, link flags.
target_link_libraries(singularity-opac::flags
INTERFACE
    MPI::MPI_CXX
    Kokkos::Kokkos
    $<${with_hdf5}:
        singularity-opac::hdf5
    >
)

if (SINGULARITY_USE_HDF5)
  target_compile_definitions(singularity-opac::flags INTERFACE
                             SPINER_USE_HDF)
  if(SINGULARITY_HDF5_INSTALL_DIR)
    list(APPEND CMAKE_PREFIX_PATH "${SINGULARITY_HDF5_INSTALL_DIR}")
  endif()
  if(HDF5_FOUND)
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