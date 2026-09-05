# Check if CMAKE_CUDA_ARCHITECTURES is set
if( NOT DEFINED CMAKE_CUDA_ARCHITECTURES )
  message( FATAL_ERROR "CMAKE_CUDA_ARCHITECTURES Must Be Set" )
endif()

# Check that only CUDA CC 6.0+ is enabled
foreach( cuda_arch ${CMAKE_CUDA_ARCHITECTURES} )
  if( cuda_arch LESS 60 )
    message(FATAL_ERROR "GauXC Requires CUDA CC 6.0+ For FP64 Atomics")
  endif()
endforeach()



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

if( GAUXC_ENABLE_NCCL )

  message( STATUS "NCCL Has Been Enabled" )
  find_package( NCCL REQUIRED )
  target_link_libraries( gauxc PUBLIC NCCL::nccl )

endif()

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

if( TARGET gauxc_cub ) # Handle the case when CUB is implicit
  target_link_libraries( gauxc PRIVATE $<BUILD_INTERFACE:gauxc_cub> )
endif()
