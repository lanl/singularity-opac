#------------------------------------------------------------------------------#
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
#------------------------------------------------------------------------------#

# Find the native EOSPAC headers and libraries.
#
#  EOSPAC_INCLUDE_DIRS - where to find eos_Interface.h, etc.
#  EOSPAC_LIBRARIES    - List of libraries when using eospac6.
#  EOSPAC_FOUND        - True if eospac found.

#TODO: Add EOSPAC_MODULES (possibly part of include dirs) and clarify which interface we need

# Look for the header file.
FIND_PATH(GANDOLF_INC_DIR NAMES gandolf.h)

# Look for the library.
FIND_LIBRARY(GANDOLF_LIB_DIR NAMES libgandolf.a)
#message("Gandolf lib dir; ${GANDOLF_LIB_DIR}")

# handle the QUIETLY and REQUIRED arguments and set EOSPAC_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(gandolf DEFAULT_MSG GANDOLF_LIB_DIR GANDOLF_INC_DIR)

# Copy the results to the output variables.
SET(GANDOLF_LIBRARY ${GANDOLF_LIB_DIR})
SET(GANDOLF_INCLUDE_DIRECTORY ${GANDOLF_INC_DIR})

MARK_AS_ADVANCED(GANDOLF_INCLUDE_DIRECTORY GANDOLF_LIBRARY)

if(GANDOLF_INCLUDE_DIRECTORY AND GANDOLF_LIBRARY)
	add_library(gandolf::gandolf STATIC IMPORTED)
	set_target_properties(gandolf::gandolf PROPERTIES
		IMPORTED_LOCATION "${GANDOLF_LIBRARY}"
		INTERFACE_INCLUDE_DIRECTORIES "${GANDOLF_INCLUDE_DIRECTORY}")
endif()

