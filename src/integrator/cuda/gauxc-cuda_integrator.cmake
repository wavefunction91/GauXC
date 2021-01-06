if( NOT TARGET CUDA::cublas )
  find_package( CUDAToolkit REQUIRED )
endif()
include( gauxc-cub )



target_sources( gauxc PRIVATE cuda/collocation_device.cu
                              cuda/xc_cuda_data.cxx
                              cuda/cuda_driver_replicated_density_incore.cxx
                              cuda/cuda_driver_replicated_density_shellbatched.cxx
                              cuda/cuda_weights.cu
                              cuda/cuda_pack_density.cu
                              cuda/cuda_eval_denvars.cu
                              cuda/cublas_extensions.cu
                              cuda/cuda_zmat.cu
                              cuda/cuda_inc_potential.cu
			      cuda/cuda_device_properties.cxx
)

target_compile_features( gauxc PRIVATE cuda_std_14 )
#target_include_directories( gauxc
#  PRIVATE
#    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/integrator/cuda>
#    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/integrator/cuda/collocation>
#)

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

target_link_libraries( gauxc PUBLIC CUDA::cublas )
target_link_libraries( gauxc PRIVATE $<BUILD_INTERFACE:gauxc_cub> )
