if( NOT TARGET CUDA::cublas )
  find_package( CUDAToolkit REQUIRED )
endif()
include( gauxc-cub )

target_sources( gauxc PRIVATE 
  # Common CUDA Utilities
  device/cuda/collocation_device.cu
  device/cuda/xc_cuda_data.cxx
  device/cuda/cuda_weights.cu
  device/cuda/cuda_pack_density.cu
  device/cuda/cuda_eval_denvars.cu
  device/cuda/cublas_extensions.cu
  device/cuda/cuda_inc_potential.cu
  device/cuda/cuda_device_properties.cxx

  # XC Specific
  device/cuda/cuda_zmat.cu

  # Drivers
  device/cuda/local_work_replicated_incore_exc_vxc.cxx

)

target_compile_features( gauxc PRIVATE cuda_std_14 )
target_compile_options( gauxc
  PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>: -Xcudafe --diag_suppress=partial_override -Xptxas -v > 
)


if( GAUXC_ENABLE_MAGMA )

  message( STATUS "MAGMA Has Been Enabled" )
  find_package( MAGMA REQUIRED )
  target_link_libraries( gauxc PUBLIC MAGMA::magma )

else()

  message( STATUS "MAGMA Has Been Explicitly Disabled" )

endif()

if(NOT GAUXC_LINK_CUDA_STATIC)
  target_link_libraries( gauxc PUBLIC CUDA::cublas )
else()
  target_link_libraries( gauxc PUBLIC CUDA::cublas_static )
endif()
target_link_libraries( gauxc PRIVATE $<BUILD_INTERFACE:gauxc_cub> )
