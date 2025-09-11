set(AUTODIFF_REQUIRED_VERSION 1.1.2)
cmake_policy(SET CMP0135 NEW)

# list(APPEND CMAKE_PREFIX_PATH "${ASTRO_THIRD_PARTY_DIR}")
find_package(
  autodiff
  ${AUTODIFF_REQUIRED_VERSION}
  NO_MODULE
  QUIET
)

if(NOT TARGET autodiff::autodiff)
  message(STATUS "Astro: Did not find AutoDiff ${AUTODIFF_REQUIRED_VERSION} installed, downloading to "
    "${ASTRO_THIRD_PARTY_DIR}")
  include(FetchContent)

  set(FETCHCONTENT_BASE_DIR "${ASTRO_THIRD_PARTY_DIR}")
  fetchcontent_declare(
    autodiff
    URL "https://github.com/autodiff/autodiff/archive/refs/tags/v${AUTODIFF_REQUIRED_VERSION}.tar.gz"
  )

  option(AUTODIFF_BUILD_TESTS OFF)
  option(AUTODIFF_BUILD_PYTHON OFF)
  option(AUTODIFF_BUILD_EXAMPLES OFF)
  option(AUTODIFF_BUILD_DOCS OFF)
  option(AUTODIFF_BUILD_PYTHON OFF)

  FetchContent_MakeAvailable(autodiff)
else()
  get_target_property(AUTODIFF_INCLUDE_DIRS
    autodiff::autodiff
    INTERFACE_INCLUDE_DIRECTORIES
  )
  message(STATUS "Astro: Found AutoDiff installed in ${AUTODIFF_INCLUDE_DIRS}")
endif()

target_compile_options(autodiff::autodiff INTERFACE
  -Wno-sign-compare
)
