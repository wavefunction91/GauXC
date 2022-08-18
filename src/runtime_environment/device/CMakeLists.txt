target_sources( gauxc PRIVATE device_runtime_environment.cxx )
if(GAUXC_ENABLE_CUDA)
  add_subdirectory( cuda )
endif()

if(GAUXC_ENABLE_HIP)
  add_subdirectory( hip )
endif()

