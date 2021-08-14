add_subdirectory(host)
target_sources( gauxc PRIVATE factory.cxx )
target_include_directories( gauxc
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
)
