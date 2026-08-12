add_executable(test_converted_data_file

    test_converted_data_file.cpp
)

target_link_libraries(test_converted_data_file
  PRIVATE
    singularity-opac::singularity-opac
    PortsofCall::PortsofCall
)

# copy converted ipcress file from analyticOpacities.ipcress
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/test_file.h5
    ${CMAKE_CURRENT_BINARY_DIR}/test_file.h5
    COPYONLY
)

add_test(
    NAME test_converted_data_file
    COMMAND test_converted_data_file
)
