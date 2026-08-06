#------------------------------------------------------------------------------#
# © 2026. Triad National Security, LLC. All rights reserved.  This
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

# ARL NOTE: Based off FindEOSPAC.cmake from singularity-eos

# Find the native GANDOLF headers and libraries.
#
#  GANDOLF_INCLUDE_DIRECTORY - Where to find gandolf.h, etc.
#  GANDOLF_LIBRARY           - Gandolf library
#  GANDOLF_FOUND             - True if gandolf found.
#
#  Current modules set GANDOLF_INC_DIR and GANDOLF_LIB_DIR and prepend the top
#  level of the Gandolf build to CMAKE_PREFIX_PATH


# if environment variables are set, use them as hints to FIND calls
set(GANDOLF_INC_DIR_HINTS "")
if(DEFINED ENV{GANDOLF_INC_DIR})
    list(APPEND GANDOLF_INC_DIR_HINTS "$ENV{GANDOLF_INC_DIR}")
endif()

set(GANDOLF_LIB_DIR_HINTS "")
if(DEFINED ENV{GANDOLF_LIB_DIR})
    list(APPEND GANDOLF_LIB_DIR_HINTS "$ENV{GANDOLF_LIB_DIR}")
endif()


# Look for the header file.
FIND_PATH(GANDOLF_INC_DIR NAMES gandolf.h HINTS ${GANDOLF_INC_DIR_HINTS})

# Look for the library.
FIND_LIBRARY(GANDOLF_LIB_DIR NAMES libgandolf.a HINTS ${GANDOLF_LIB_DIR_HINTS})

# handle the QUIETLY and REQUIRED arguments and set GANDOLF_FOUND to TRUE if
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
