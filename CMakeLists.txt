# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                        #
#                                                                                                 #
# The Astro project is distributed under the BSD 2-Clause License.                                #
#                                                                                                 #
# Davide Stocco                                                                 Enrico Bertolazzi #
# University of Trento                                                       University of Trento #
# e-mail: davide.stocco@unitn.it                               e-mail: enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

cmake_minimum_required(VERSION 3.14)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(AstroVersion)
get_version_from_git()

project(Astro
  VERSION ${ASTRO_VERSION}
  DESCRIPTION "A C++ library for space dynamics"
  HOMEPAGE_URL "https://github.com/StoccoDavide/Astro"
  LANGUAGES CXX
)

if(NOT ASTRO_VERSION_PRERELEASE)
  set(ASTRO_FULL_VERSION "${ASTRO_VERSION}")
  write_version_file()
else()
  set(ASTRO_FULL_VERSION "${ASTRO_VERSION}-${ASTRO_VERSION_PRERELEASE}")
endif()
message(STATUS "Astro: version ${ASTRO_FULL_VERSION}")

set(ASTRO_THIRD_PARTY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third_party")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

if(NOT CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/lib")
  message(STATUS  "Astro: CMAKE_INSTALL_PREFIX not specified, defaulting to ${CMAKE_CURRENT_SOURCE_DIR}")
endif()

# Dependencies
include(AstroEigen3)
include(AstroSandals)
include(AstroAutoDiff)

# Library definition
add_library(Astro INTERFACE)
add_library(Astro::Astro ALIAS Astro)

target_link_libraries(Astro INTERFACE Eigen3::Eigen Sandals::Sandals autodiff::autodiff)

target_include_directories(Astro INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

# INFO: Ensure this options are only available to developers. In theory we could leave them available,
# but the call to `target_sources` breaks `fetchcontent` compatibility due to absolute paths being
# added to `INTERFACE_SOURCES`. I tried solving it, but it seems to be poorly documented, supported, etc.
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
    message(STATUS
      "Astro: "
      "CMAKE_BUILD_TYPE not specified, defaulting to ${CMAKE_BUILD_TYPE}"
    )
  endif()

  option(BUILD_TESTS "Build tests" OFF)
  option(BUILD_BENCHMARKS "Build benchmarks" OFF)
  option(BUILD_EXAMPLES "Build examples" OFF)
  option(ASTRO_ENABLE_PLOTTING "Enable plotting" OFF)

  if(ASTRO_ENABLE_PLOTTING)
    include(AstroRoot)
    target_link_libraries(Astro INTERFACE ROOT::Core ROOT::Gpad ROOT::Geom ROOT::Graf3d)
  endif()

  file(GLOB_RECURSE HEADER_FILES_HH "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hh")
  file(GLOB_RECURSE HEADER_FILES_HXX "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hxx")
  foreach(HEADER_FILE IN LISTS HEADER_FILES_HH HEADER_FILES_HXX)
    file(
      RELATIVE_PATH
      RELATIVE_PATH_HEADER_FILE
      ${CMAKE_CURRENT_SOURCE_DIR} ${HEADER_FILE}
    )
    list(APPEND RELATIVE_HEADER_FILES ${RELATIVE_PATH_HEADER_FILE})
  endforeach()
  target_sources(Astro INTERFACE $<BUILD_INTERFACE:${RELATIVE_HEADER_FILES}>)

  if(BUILD_TESTS)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tests)
  endif()

  if(BUILD_BENCHMARKS)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/benchmarks)
  endif()

  if(BUILD_EXAMPLES)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/examples)
  endif()

  if(ASTRO_ENABLE_PLOTTING)
    target_compile_definitions(Astro INTERFACE ASTRO_ENABLE_PLOTTING)
  endif()
endif()

if(MSVC)
  target_compile_options(Astro INTERFACE /W4 /WX)
else()
  target_compile_options(Astro INTERFACE -Wall -Wextra -pedantic -Werror)
endif()

# Installation
install(TARGETS Astro EXPORT AstroTargets)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION include)

include(CMakePackageConfigHelpers)

export(EXPORT AstroTargets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/AstroTargets.cmake
  NAMESPACE Astro::
)

install(EXPORT AstroTargets
  FILE AstroTargets.cmake
  NAMESPACE Astro::
  DESTINATION lib/cmake/Astro
)

write_basic_package_version_file(
  AstroConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/AstroConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/AstroConfig.cmake
  INSTALL_DESTINATION lib/cmake/Astro
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/AstroConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/AstroConfigVersion.cmake
  DESTINATION lib/cmake/Astro
)
