cmake_minimum_required(VERSION 3.15)
project(exotic-quadrature-experiments)

# Download and unpack TASMANIAN.
include(FetchContent)
set(Tasmanian_ENABLE_PYTHON ON CACHE INTERNAL "Turn on TASMANIAN Python wrappers.")
FetchContent_Declare(
  Tasmanian
  GIT_REPOSITORY    https://github.com/ORNL/TASMANIAN.git
  GIT_TAG           master
)
FetchContent_MakeAvailable(tasmanian)

# Build Dependencies.
add_subdirectory("${CMAKE_SOURCE_DIR}/src")
add_subdirectory("${CMAKE_SOURCE_DIR}/scripts")
