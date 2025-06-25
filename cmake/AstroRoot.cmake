if (ASTRO_ENABLE_PLOTTING)
  set(ROOT_REQUIRED_VERSION 6.34.08)
  cmake_policy(SET CMP0135 NEW)

  # list(APPEND CMAKE_PREFIX_PATH "${ASTRO_THIRD_PARTY_DIR}")
  set(ROOT_REQUIRED_COMPONENTS Core Gpad)
  find_package(
    ROOT
    ${ROOT_REQUIRED_VERSION}
    REQUIRED COMPONENTS ${ROOT_REQUIRED_COMPONENTS}
    NO_MODULE
    QUIET
  )

  get_target_property(ROOT_INCLUDE_DIRS
    ROOT::Core
    INTERFACE_INCLUDE_DIRECTORIES
  )

  if(NOT TARGET ROOT::Core OR NOT TARGET ROOT::Gpad)
    message(STATUS "Astro: Did not find Root ${ROOT_REQUIRED_VERSION} installed, please "
      "install it manually via your package manager.")
  else()
    message(STATUS "Astro: Found Root installed in ${ROOT_INCLUDE_DIRS}")
  endif()
else()
  add_library(ROOT::Core INTERFACE IMPORTED)
  add_library(ROOT::Gpad INTERFACE IMPORTED)
  message(STATUS "Astro: Root plotting disabled")
endif()
