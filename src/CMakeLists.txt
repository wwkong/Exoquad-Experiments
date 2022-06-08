set(Exoquad_source_files ExpGauss.jl
                         ExoquadUtils.py)

# Copies each source file to the build directory
foreach(Exoquad_source_file ${Exoquad_source_files})
  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${Exoquad_source_file}"
    COMMAND "${CMAKE_COMMAND}"
    ARGS    -E copy "${CMAKE_CURRENT_SOURCE_DIR}/${Exoquad_source_file}" "${CMAKE_CURRENT_BINARY_DIR}/${Exoquad_source_file}"
    COMMENT "Copying ${CMAKE_CURRENT_SOURCE_DIR}/${Exoquad_source_file}")
  list(APPEND Exoquad_pre_source_files "${CMAKE_CURRENT_BINARY_DIR}/${Exoquad_source_file}")
endforeach()

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/../Project.toml"
  COMMAND "${CMAKE_COMMAND}"
  ARGS    -E copy "${CMAKE_CURRENT_SOURCE_DIR}/Project.toml" "${CMAKE_CURRENT_BINARY_DIR}/../Project.toml"
  COMMENT "Copying ${CMAKE_CURRENT_SOURCE_DIR}/Project.toml")
list(APPEND Exoquad_pre_source_files "${CMAKE_CURRENT_BINARY_DIR}/../Project.toml")

add_custom_target(Exoquad_source ALL DEPENDS "${Exoquad_pre_source_files}")
