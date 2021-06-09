


# Easier generator expressions
set(build_debug "$<CONFIG:Debug>")
set(build_release "$<CONFIG:Release>")
set(cxx_lang "$<COMPILE_LANGUAGE:CXX>")
set(cxx_xl "$<COMPILE_LANG_AND_ID:CXX,XL>")
set(with_hdf5 "$<BOOL:${SINGULARITY_USE_HDF5}>")
set(with_kokkos "$<BOOL:${SINGULARITY_USE_KOKKOS}>")
set(with_cuda "$<BOOL:${SINGULARITY_USE_CUDA}>")
set(hide_more_warn "$<BOOL:${SINGULARITY_HIDE_MORE_WARNINGS}>")
set(better_debug "$<BOOL:${SINGULARITY_BETTER_DEBUG_FLAGS}>")

# xl fix
target_compile_options(${PROJECT_NAME} INTERFACE
                        $<${cxx_xl}:
                            "-std=c++1y;-qxflag=disable__cplusplusOverride"
                        >
)
target_link_options(${PROJECT_NAME} INTERFACE
                        $<${cxx_xl}:
                            "-std=c++1y;-qxflag=disable__cplusplusOverride"
                        >
)

# Base Include directories
target_include_directories(${PROJECT_NAME}
  INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/utils/herumi-fmath>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/utils>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>

)

target_compile_options(${PROJECT_NAME}
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
target_link_libraries(${PROJECT_NAME}
INTERFACE
    MPI::MPI_CXX
    $<${with_kokkos}:
        Kokkos::Kokkos
    >
    $<${with_hdf5}:
        ${PROJECT_NAME}::hdf5
    >
)
