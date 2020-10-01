
target_sources( gauxc PRIVATE hip/collocation_device.cxx
                              hip/xc_hip_data.cxx
                              hip/xc_hip_util.cxx
                              hip/hip_weights.cxx
                              hip/hip_pack_density.cxx
                              hip/hip_eval_denvars.cxx
                              hip/hipblas_extensions.cxx
                              hip/hip_zmat.cxx
                              hip/hip_inc_potential.cxx
)

target_include_directories( gauxc
  PRIVATE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/integrator/hip>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/integrator/hip/collocation>
)

#find_package( hipblas REQUIRED )
#target_link_libraries( gauxc PUBLIC roc::hipblas )

add_library( hipblas INTERFACE IMPORTED )
set_target_properties( hipblas PROPERTIES
  INTERFACE_LINK_LIBRARIES "-L/opt/rocm-3.8.0/hipblas/lib -lhipblas"
  INTERFACE_INCLUDE_DIRECTORIES "/opt/rocm-3.8.0/hipblas/include"
)

target_link_libraries( gauxc PUBLIC hipblas )
