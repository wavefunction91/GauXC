target_sources( gauxc PRIVATE runtime_environment.cxx )
target_include_directories( gauxc
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
)
if(GAUXC_ENABLE_DEVICE)
  add_subdirectory(device)
endif()
