target_sources( gauxc PRIVATE 
  device_reduction_driver.cxx
)

if( GAUXC_ENABLE_NCCL )
  target_sources( gauxc PRIVATE nccl_reduction_driver.cxx )
  message( STATUS "NCCL Has Been Enabled" )
  find_package( NCCL REQUIRED )
  target_link_libraries( gauxc PUBLIC NCCL::nccl )
endif()

