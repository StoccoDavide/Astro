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

file(GLOB_RECURSE ORBIT_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/Orbit/test_*.cc")
add_executable(test_orbit ${ORBIT_TESTS})
target_link_libraries(test_orbit PRIVATE Astro)

file(GLOB_RECURSE ORBIT_PLOTTING "${CMAKE_CURRENT_SOURCE_DIR}/Plotting/test_*.cc")
add_executable(test_plotting ${ORBIT_PLOTTING})
target_link_libraries(test_plotting PRIVATE Astro)# ROOT::Core ROOT::Gpad ROOT::Geom)
