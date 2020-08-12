target_sources( gauxc
  PRIVATE
    integrator_common.cxx
    host/xc_host_util.cxx
    host/host_weights.cxx
    host/host_collocation.cxx
    host/host_zmat.cxx
    host/blas.cxx
)

target_include_directories( gauxc
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/integrator>
)

if( GAUXC_ENABLE_DEVICE )

  if( GAUXC_ENABLE_CUDA )
    include( cuda/gauxc-cuda_integrator.cmake )
  endif()

  get_target_property( DEVICE_TARGET_DEPENDS gauxc_device_integrator
    INTERFACE_LINK_LIBRARIES ) 

  target_link_libraries( gauxc PUBLIC 
    $<TARGET_OBJECTS:gauxc_device_integrator> 
    ${DEVICE_TARGET_DEPENDS}
  )

  add_dependencies( gauxc gauxc_device_integrator )

endif()
