target_sources( gauxc PRIVATE 
  # Drivers
  device/local_work_replicated_shellbatched_exc_vxc.cxx

  # Interfaces
  device/incore_xc_device_integrator.cxx
  device/shellbatched_xc_device_integrator.cxx
)

if( GAUXC_ENABLE_CUDA )
  include( device/cuda/gauxc-cuda.cmake )
endif()


if( GAUXC_ENABLE_HIP )
  include( device/hip/gauxc-hip.cmake )
endif()
