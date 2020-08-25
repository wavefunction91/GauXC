find_package( CUDAToolkit REQUIRED )
find_package( MAGMA       REQUIRED )


#add_library( gauxc_device_integrator STATIC
#  #cuda/collocation_kernels.cu
#  cuda/collocation_device.cu
#  cuda/xc_cuda_data.cxx
#  cuda/xc_cuda_util.cxx
#)
#target_compile_features( gauxc_device_integrator PUBLIC cuda_std_14 cxx_std_14 )
#target_link_libraries( gauxc_device_integrator PUBLIC MAGMA::magma CUDA::cudart 
#                       CUDA::cublas gauxc_integrator_common )
#
#target_include_directories( gauxc_device_integrator
#  PRIVATE
#    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
#    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
#    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
#    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/integrator>
#    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/integrator/cuda>
#    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/integrator/cuda/collocation>
#)
#
#target_compile_options( gauxc_device_integrator 
#  PRIVATE
#    $<$<COMPILE_LANGUAGE:CXX>: -Wall -Wextra -Wpedantic -Wnon-virtual-dtor>
#    $<$<COMPILE_LANGUAGE:CUDA>: -Xcudafe --diag_suppress=partial_override -Xptxas -v > 
#)
#


target_sources( gauxc PRIVATE cuda/collocation_device.cu
                              cuda/xc_cuda_data.cxx
                              cuda/xc_cuda_util.cxx
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


