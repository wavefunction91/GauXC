target_sources( gauxc PRIVATE 
  local_device_work_driver.cxx
  local_device_work_driver_pimpl.cxx

  scheme1_base.cxx
  scheme1_data_base.cxx
)

if( GAUXC_ENABLE_MAGMA )
  find_package( MAGMA REQUIRED )
  target_link_libraries( gauxc PUBLIC MAGMA::magma )

  target_sources( gauxc PRIVATE 
    scheme1_magma_base.cxx
    scheme1_magma_data_base.cxx
  )
endif()

if(GAUXC_ENABLE_CUDA)

  add_subdirectory( cuda )
endif()

if(GAUXC_ENABLE_HIP)
  add_subdirectory( hip )
endif()
