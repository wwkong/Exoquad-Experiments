message(STATUS "Subproject: TASMANIAN...")
# Download and unpack gflags at configure time
configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/tasmanian.CMakeLists.txt
	${CMAKE_BINARY_DIR}/tasmanian-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tasmanian-download)
if(result)
  message(FATAL_ERROR "CMake step for TASMANIAN failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tasmanian-download)
if(result)
  message(FATAL_ERROR "Build step for TASMANIAN failed: ${result}")
endif()

add_subdirectory(
	${CMAKE_BINARY_DIR}/tasmanian-src
	${CMAKE_BINARY_DIR}/tasmanian-build)
message(STATUS "Subproject: TASMANIAN...DONE")
