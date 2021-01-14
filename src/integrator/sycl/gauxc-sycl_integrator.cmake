include( gauxc-eigen3 )
target_sources( gauxc PRIVATE sycl/collocation_device.cpp
                              sycl/xc_sycl_data.cxx
                              sycl/xc_sycl_util.cxx
                              sycl/sycl_weights.cpp
                              sycl/sycl_pack_density.cpp
                              sycl/sycl_eval_denvars.cpp
                              sycl/mklsycl_extensions.cpp
                              sycl/sycl_zmat.cpp
                              sycl/sycl_inc_potential.cpp
)

target_include_directories( gauxc
  PRIVATE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/integrator/sycl>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/integrator/sycl/collocation>
)

target_link_libraries( gauxc PRIVATE $<BUILD_INTERFACE:Eigen3::Eigen> )

# Handle the case when HOST is disabled but we still need oneMKL
if( NOT GAUXC_ENABLE_HOST )
  find_package( LAPACK REQUIRED COMPONENTS sycl )
  target_link_libraries( gauxc PUBLIC LAPACK::LAPACK )
endif()
