#
# GauXC Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of
# any required approvals from the U.S. Dept. of Energy). All rights reserved.
#
# See LICENSE.txt for details
#
if( NOT TARGET CUDA::cublas )
  find_package( CUDAToolkit REQUIRED )
endif()

target_compile_features( gauxc PRIVATE cuda_std_14 )
target_sources( gauxc PRIVATE cuda_backend.cxx )

if(BUILD_SHARED_LIBS)
  target_link_libraries( gauxc PUBLIC CUDA::cudart CUDA::cublas )
else()
  target_link_libraries( gauxc PUBLIC CUDA::cudart_static CUDA::cublas_static )
endif()
