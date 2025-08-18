# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                        #
#                                                                                                 #
# The Astro project is distributed under the BSD 2-Clause License.                                #
#                                                                                                 #
# Davide Stocco                                                                 Enrico Bertolazzi #
# University of Trento                                                       University of Trento #
# e-mail: davide.stocco@unitn.it                               e-mail: enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

include(AstroCatch2)

file(GLOB_RECURSE ANOMALY_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/Anomaly/test_*.cc")
add_executable(test_anomaly ${ANOMALY_TESTS})
target_link_libraries(test_anomaly PRIVATE Astro Catch2::Catch2WithMain)

if(ASTRO_ENABLE_PLOTTING)
  #file(GLOB_RECURSE ORBIT_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/Orbit/test_*.cc")
  #add_executable(test_orbit ${ORBIT_TESTS})
  #target_link_libraries(test_orbit PRIVATE Astro)
#
  #file(GLOB_RECURSE ORBIT_PLOTTING "${CMAKE_CURRENT_SOURCE_DIR}/Plotting/test_*.cc")
  #add_executable(test_plotting ${ORBIT_PLOTTING})
  #target_link_libraries(test_plotting PRIVATE Astro)
#
  #file(GLOB_RECURSE ORBIT_BODY "${CMAKE_CURRENT_SOURCE_DIR}/Body/test_*.cc")
  #add_executable(test_body ${ORBIT_BODY})
  #target_link_libraries(test_body PRIVATE Astro)
#
  #file(GLOB_RECURSE ORBIT_MAGNETICFIELD "${CMAKE_CURRENT_SOURCE_DIR}/MagneticField/test_*.cc")
  #add_executable(test_magneticfield ${ORBIT_MAGNETICFIELD})
  #target_link_libraries(test_magneticfield PRIVATE Astro)

  file(GLOB_RECURSE TEST_SUNEARTHMOON "${CMAKE_CURRENT_SOURCE_DIR}/SunEarthMoon/test_*.cc")
  add_executable(test_sunearthmoon ${TEST_SUNEARTHMOON})
  target_link_libraries(test_sunearthmoon PRIVATE Astro)
else()
  message(STATUS "Astro: Root plotting disabled, skipping tests involving plotting")
endif()
