include( gauxc-gau2grid     )
target_sources( gauxc PRIVATE host/xc_host_util.cxx
                              host/host_weights.cxx
                              host/host_collocation.cxx
                              host/host_zmat.cxx
                              host/blas.cxx
)

set( GAUXC_LAPACK_OPTIONAL_COMPONENTS )
set( GAUXC_LAPACK_REQUIRED_COMPONENTS )

if( GAUXC_ENABLE_SYCL )
  set( GAUXC_LAPACK_REQUIRED_COMPONENTS sycl )
endif()

find_package( LAPACK REQUIRED COMPONENTS ${GAUXC_LAPACK_REQUIRED_COMPONENTS} )


target_link_libraries( gauxc PUBLIC LAPACK::LAPACK )

if( GAUXC_ENABLE_GAU2GRID )
  target_link_libraries( gauxc PUBLIC gau2grid::gg )
endif()

