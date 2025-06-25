set(SANDALS_REQUIRED_VERSION 0.0.0)
cmake_policy(SET CMP0135 NEW)

# list(APPEND CMAKE_PREFIX_PATH "${ASTRO_THIRD_PARTY_DIR}")
find_package(
  Sandals
  ${SANDALS_REQUIRED_VERSION}
  NO_MODULE
  QUIET
)

if(NOT TARGET Sandals::Sandals)
  message(STATUS "Astro: Did not find Sandals ${SANDALS_REQUIRED_VERSION} installed, downloading to "
    "${ASTRO_THIRD_PARTY_DIR}")
  include(FetchContent)

  set(FETCHCONTENT_BASE_DIR "${ASTRO_THIRD_PARTY_DIR}")
  fetchcontent_declare(
    Sandals
    GIT_REPOSITORY "https://github.com/StoccoDavide/Sandals"
    GIT_SHALLOW    TRUE
    GIT_TAG        main
  )

  option(SANDALS_BUILD_TESTS OFF)
  option(SANDALS_BUILD_EXAMPLES OFF)
  option(SANDALS_BUILD_BENCHMARKS OFF)
  option(SANDALS_ENABLE_PLOTTING OFF)

  fetchcontent_makeavailable(Sandals)
else()
  get_target_property(SANDALS_INCLUDE_DIRS
    Sandals::Sandals
    INTERFACE_INCLUDE_DIRECTORIES
  )
  message(STATUS "Astro: Found Sandals installed in ${SANDALS_INCLUDE_DIRS}")
endif()
