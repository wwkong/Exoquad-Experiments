# Initialization.
cmake_minimum_required(VERSION 3.0)
project(exotic-quadrature-experiments)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Build TASMANIAN.
execute_process(
  COMMAND ${CMAKE_COMMAND} ${CMAKE_SOURCE_DIR}/cmake
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/TASMANIAN
)
execute_process(
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR}/external/TASMANIAN
)

# # Build Dependencies.
# add_subdirectory(cmake)
# add_subdirectory(julia)
# add_subdirectory(python)

