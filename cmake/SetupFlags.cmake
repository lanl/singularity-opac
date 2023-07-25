# Easier generator expressions
set(build_debug "$<CONFIG:Debug>")
set(build_release "$<CONFIG:Release>")
set(cxx_lang "$<COMPILE_LANGUAGE:CXX>")
set(cxx_xl "$<COMPILE_LANG_AND_ID:CXX,XL>")
set(with_fmath "$<BOOL:${SINGULARITY_USE_FMATH}>")
set(with_hdf5 "$<BOOL:${SINGULARITY_USE_HDF5}>")
set(with_mpi "$<BOOL:${SINGULARITY_USE_MPI}>")
set(with_kokkos "$<BOOL:${SINGULARITY_USE_KOKKOS}>")
set(without_kokkos "$<NOT:${with_kokkos}>")
set(with_cuda "$<BOOL:${SINGULARITY_USE_CUDA}>")
set(without_cuda "$<NOT:${with_cuda}>")
set(hide_more_warn "$<BOOL:${SINGULARITY_HIDE_MORE_WARNINGS}>")
set(better_debug "$<BOOL:${SINGULARITY_BETTER_DEBUG_FLAGS}>")

# xl fix
target_compile_options(
  singularity-opac INTERFACE $<${cxx_xl}:
                             "-std=c++1y;-qxflag=disable__cplusplusOverride" >)
target_link_options(singularity-opac INTERFACE $<${cxx_xl}:
                    "-std=c++1y;-qxflag=disable__cplusplusOverride" >)

# Base Include directories target_include_directories( singularity-opac
# INTERFACE $<${without_kokkos}:
# $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/utils/herumi-fmath> >
# $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/utils>
# $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>)

target_compile_options(
  singularity-opac
  INTERFACE $<${with_kokkos}:
            $<${with_cuda}:
            $<${cxx_lang}:
            "--expt-relaxed-constexpr"
            $<${hide_more_warn}:
            "-Xcudafe;--diag_suppress=esa_on_defaulted_function_ignored"
            > # hide_more_warn
            > # cxx_lang
            $<${build_release}:
            "-use_fast_math"
            > # build_release
            $<${build_debug}:
            $<${better_debug}:
            $<${cxx_lang}:
            "-G;-lineinfo"
            > # cxx_lang
            > # better_debug
            > # build_debug
            > # with_cuda
            > # with_kokkos
)

target_compile_definitions(
  singularity-opac
  INTERFACE $<${with_kokkos}:
            PORTABILITY_STRATEGY_KOKKOS
            >
            $<${without_kokkos}:
            $<${with_fmath}:
            SINGULARITY_USE_FMATH
            >
            >)
