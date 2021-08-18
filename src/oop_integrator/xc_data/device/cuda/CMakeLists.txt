if( NOT TARGET CUDA::cublas )
  find_package( CUDAToolkit REQUIRED )
endif()

target_compile_features( gauxc PRIVATE cuda_std_14 )
target_sources( gauxc PRIVATE xc_cuda_aos_data_base.cxx )
target_link_libraries( gauxc PUBLIC CUDA::cudart CUDA::cublas )
