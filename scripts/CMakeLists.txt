# Main Julia Scripts
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/../res")
set(Exoquad_res_folder "\"${CMAKE_CURRENT_BINARY_DIR}/../res\"")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/beylkin_antenna.in.jl" "${CMAKE_CURRENT_BINARY_DIR}/beylkin_antenna.jl")

# Key variables
set(Tasmanian_python_import "sys.path.append(\"${Tasmanian_BINARY_DIR}/Python\")")
set(Exoquad_src_python_import "sys.path.append(\"${CMAKE_CURRENT_BINARY_DIR}/../src\")")
set(Exoquad_script_python_import "sys.path.append(\"${CMAKE_CURRENT_BINARY_DIR}/../script\")")

# Main Python Scripts
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GenerateData.in.py" "${CMAKE_CURRENT_BINARY_DIR}/GenerateData.py")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GeneratePlots_1D.in.py" "${CMAKE_CURRENT_BINARY_DIR}/GeneratePlots_1D.py")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GeneratePlots_nD1.in.py" "${CMAKE_CURRENT_BINARY_DIR}/GeneratePlots_nD1.py")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GeneratePlots_nD2.in.py" "${CMAKE_CURRENT_BINARY_DIR}/GeneratePlots_nD2.py")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GeneratePlots_nD3.in.py" "${CMAKE_CURRENT_BINARY_DIR}/GeneratePlots_nD3.py")

# Auxiliary scripts
set(Exoquad_aux_script_files AntennaExp.py
                             SincExp.py
                             AbsSinExp.py
                             SqrSincExp.py
                             SumSincExp.py
                             QuadSincExp.py)
foreach(Exoquad_script_file ${Exoquad_aux_script_files})
  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${Exoquad_script_file}"
    COMMAND "${CMAKE_COMMAND}"
    ARGS    -E copy "${CMAKE_CURRENT_SOURCE_DIR}/${Exoquad_script_file}" "${CMAKE_CURRENT_BINARY_DIR}/${Exoquad_script_file}"
    COMMENT "Copying ${CMAKE_CURRENT_SOURCE_DIR}/${Exoquad_script_file}")
  list(APPEND Exoquad_pre_script_files "${CMAKE_CURRENT_BINARY_DIR}/${Exoquad_script_file}")
endforeach()

add_custom_target(Exoquad_script ALL DEPENDS "${Exoquad_pre_script_files}")
