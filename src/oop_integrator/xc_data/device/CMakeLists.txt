target_sources( gauxc PRIVATE xc_device_stack_data.cxx xc_device_aos_data.cxx )
if(GAUXC_ENABLE_CUDA)
  add_subdirectory( cuda )
endif()
