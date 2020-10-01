
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

# XXX: Need to create a good find module for this
#find_package( hipblas REQUIRED )
#target_link_libraries( gauxc PUBLIC roc::hipblas )

if( NOT HIPBLAS_LIBRARIES OR NOT HIPBLAS_INCLUDE_DIRECTORES )
  message( FATAL_ERROR "hipBLAS is not automatically discovered, you must set both HIPBLAS_LIBRARIES and HIPBLAS_INCLUDE_DIRECTORIES" )
endif()
add_library( hipblas INTERFACE IMPORTED )
set_target_properties( hipblas PROPERTIES
  INTERFACE_LINK_LIBRARIES      "${HIPBLAS_LIBRARIES}"
  INTERFACE_INCLUDE_DIRECTORIES "${HIPBLAS_INCLUDE_DIRECTORIES}"
)

target_link_libraries( gauxc PUBLIC hipblas )
