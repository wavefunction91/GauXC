find_package( hipblas REQUIRED )
#include( gauxc-cub )

target_sources( gauxc PRIVATE 
  # Common HIP Utilities
  device/hip/collocation_device.hip
  device/hip/xc_hip_data.cxx
  device/hip/hip_weights.hip
  device/hip/hip_pack_density.hip
  device/hip/hip_eval_denvars.hip
  device/hip/hipblas_extensions.hip
  device/hip/hip_inc_potential.hip
  device/hip/hip_device_properties.cxx

  # XC Specific
  device/hip/hip_zmat.hip

  # Drivers
  device/hip/local_work_replicated_incore_exc_vxc.cxx

)

#target_compile_features( gauxc PRIVATE hip_std_14 )
#target_compile_options( gauxc
#  PRIVATE
#    $<$<COMPILE_LANGUAGE:HIP>: -Xhipfe --diag_suppress=partial_override -Xptxas -v > 
#)


if( GAUXC_ENABLE_MAGMA )

  message( STATUS "MAGMA Has Been Enabled" )
  find_package( MAGMA REQUIRED )
  target_link_libraries( gauxc PUBLIC MAGMA::magma )

else()

  message( STATUS "MAGMA Has Been Explicitly Disabled" )

endif()

target_link_libraries( gauxc PUBLIC roc::hipblas )
#target_link_libraries( gauxc PRIVATE $<BUILD_INTERFACE:gauxc_cub> )
