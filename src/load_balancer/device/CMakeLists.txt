target_sources( gauxc PRIVATE 
  load_balancer_device_factory.cxx
)

if( GAUXC_ENABLE_CUDA )
  add_subdirectory( cuda )
endif()

if( GAUXC_ENABLE_HIP )
  add_subdirectory( hip )
endif()
