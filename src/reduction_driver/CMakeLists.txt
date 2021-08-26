target_sources( gauxc PRIVATE 
  reduction_driver.cxx
  reduction_driver_factory.cxx
  reduction_driver_impl.cxx
)


target_include_directories( gauxc
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
)

add_subdirectory(host)
