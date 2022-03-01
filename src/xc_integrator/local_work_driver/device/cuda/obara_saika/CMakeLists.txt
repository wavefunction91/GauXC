set( GAUXC_OBARA_SAIKA_CUDA_SRC
     src/integral_0.cu
     src/integral_1.cu
     src/integral_2.cu
     src/integral_0_0.cu
     src/integral_1_0.cu
     src/integral_1_1.cu
     src/integral_2_0.cu
     src/integral_2_1.cu
     src/integral_2_2.cu
     src/obara_saika_integrals.cu
     src/chebyshev_boys_computation.cu
)
target_sources( gauxc PRIVATE ${GAUXC_OBARA_SAIKA_CUDA_SRC} )
target_include_directories( gauxc PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
)
