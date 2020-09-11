find_package( CUDAToolkit REQUIRED )
find_package( MAGMA       REQUIRED )



target_sources( gauxc PRIVATE cuda/collocation_device.cu
                              cuda/xc_cuda_data.cxx
                              cuda/xc_cuda_util.cxx
                              cuda/cuda_weights.cu
                              cuda/cuda_pack_density.cu
                              cuda/cuda_eval_denvars.cu
                              cuda/cublas_extensions.cu
                              cuda/cuda_zmat.cu
                              cuda/cuda_inc_potential.cu
)

target_compile_features( gauxc PRIVATE cuda_std_14 )
target_include_directories( gauxc
  PRIVATE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/integrator/cuda>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/integrator/cuda/collocation>
)

target_compile_options( gauxc
  PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>: -Xcudafe --diag_suppress=partial_override -Xptxas -v > 
)

target_link_libraries( gauxc PUBLIC MAGMA::magma CUDA::cublas )


