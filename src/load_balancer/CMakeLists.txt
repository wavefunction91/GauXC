target_sources( gauxc PRIVATE 
  load_balancer.cxx 
  load_balancer_impl.cxx 
  load_balancer_factory.cxx

  host/load_balancer_host_factory.cxx
  host/replicated_host_load_balancer.cxx 
)

target_include_directories( gauxc
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
)

if( GAUXC_ENABLE_DEVICE )
  add_subdirectory( device )
endif()

