cmake_minimum_required(VERSION 3.15)
project(experiments-scripts)

# Main Script
set(Tasmanian_python_import "sys.path.append(\"${Tasmanian_BINARY_DIR}/Python\")")
set(Exoquad_python_import "sys.path.append(\"${CMAKE_CURRENT_BINARY_DIR}/../src\")")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/exoquad.in.py" "${CMAKE_CURRENT_BINARY_DIR}/exoquad.py")

# Auxiliary scripts
set(Exoquad_script_files beylkin_antenna.jl
                         exoquad_1D.py
                         exoquad_nD.py)

# Copies each source file to the build directory
foreach(Exoquad_script_file ${Exoquad_script_files})
  add_custom_command(
    PRE_BUILD
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${Exoquad_script_file}"
    COMMAND "${CMAKE_COMMAND}"
    ARGS    -E copy "${CMAKE_CURRENT_SOURCE_DIR}/${Exoquad_script_file}" "${CMAKE_CURRENT_BINARY_DIR}/${Exoquad_script_file}"
    COMMENT "Copying ${CMAKE_CURRENT_SOURCE_DIR}/${Exoquad_script_file}")
  list(APPEND Exoquad_pre_script_files "${CMAKE_CURRENT_BINARY_DIR}/${Exoquad_script_file}")
endforeach()
add_custom_target(Exoquad_script ALL DEPENDS "${Exoquad_pre_script_files}")
