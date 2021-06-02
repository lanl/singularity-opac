# Unit tests
if (SINGULARITY_BUILD_TESTS)

  include(CTest)

  # Get catch2
  message(STATUS "Fetching Catch2 as needed")
  Include(FetchContent)
  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v2.13.1)
  FetchContent_MakeAvailable(Catch2)
  list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/contrib)

  # Build tests
  message(STATUS "Configuring unit tests")
  add_subdirectory(test)
endif()