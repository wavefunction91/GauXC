add_subdirectory(integrator_util)
add_subdirectory(local_work_driver)
add_subdirectory(replicated)
add_subdirectory(xc_data)

target_include_directories( gauxc
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
)
