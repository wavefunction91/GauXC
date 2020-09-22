find_package( HIP REQUIRED )

hip_add_library( gauxc_hip  hip/collocation_device.cxx
                            hip/xc_hip_data.cxx
                            hip/xc_hip_util.cxx
                            hip/hip_weights.cxx
                            hip/hip_pack_density.cxx
                            hip/hip_eval_denvars.cxx
                            hip/hipblas_extensions.cxx
                            hip/hip_zmat.cxx
                            hip/hip_inc_potential.cxx
)

target_include_directories( gauxc_hip
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  PRIVATE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/integrator/hip>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/integrator/hip/collocation>
)
target_link_libraries( gauxc_hip PUBLIC ExchCXX::exchcxx )

#target_compile_features( gauxc PRIVATE hip_std_14 )
target_include_directories( gauxc
  PRIVATE
    ${HIP_INCLUDE_DIRECTORIES}
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/integrator/hip>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/integrator/hip/collocation>
)


#target_link_libraries( gauxc PUBLIC HIP::hipblas )
