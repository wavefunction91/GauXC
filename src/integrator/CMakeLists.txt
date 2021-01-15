# Common Integrator Utilities
target_sources( gauxc PRIVATE integrator_common.cxx )
target_include_directories( gauxc
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/integrator>
)

# Host Integrator Utilities
if( GAUXC_ENABLE_HOST )
  include( host/gauxc-host_integrator.cmake )
endif()

if( GAUXC_ENABLE_CUDA )
  include( cuda/gauxc-cuda_integrator.cmake )
endif()
