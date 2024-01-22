#=======================================
# Compiler logic
#=======================================

# set C++ options
enable_language(CXX)
include(CMakeDetermineCXXCompiler)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Fortran options
if(SINGULARITY_USE_FORTRAN)
  enable_language(Fortran)
  include(CMakeDetermineFortranCompiler) #NOTE: this can be annoying with Cray machines
endif()

