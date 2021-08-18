target_sources( gauxc PRIVATE 
  local_device_work_driver.cxx
  local_device_work_driver_pimpl.cxx
)

if(GAUXC_ENABLE_CUDA)
  add_subdirectory( cuda )
endif()
