find_package( ExchCXX QUIET )
if( NOT ${ExchCXX_FOUND} )

  include( gauxc-dep-versions )

  message( STATUS "Could not find ExchCXX... Building" )
  message( STATUS "EXCHCXX REPO = ${GAUXC_EXCHCXX_REPOSITORY}" )
  message( STATUS "EXCHCXX REV  = ${GAUXC_EXCHCXX_REVISION}"   )

  set( EXCHCXX_ENABLE_CUDA  ${GAUXC_ENABLE_CUDA} CACHE BOOL "" )
  set( EXCHCXX_ENABLE_SYCL  ${GAUXC_ENABLE_SYCL} CACHE BOOL "" )
  set( EXCHCXX_ENABLE_TESTS OFF                  CACHE BOOL "" )

  FetchContent_Declare(
    exchcxx
    GIT_REPOSITORY ${GAUXC_EXCHCXX_REPOSITORY} 
    GIT_TAG        ${GAUXC_EXCHCXX_REVISION} 
  )

  FetchContent_MakeAvailable( exchcxx )


else()

  if( ${GAUXC_ENABLE_CUDA} AND NOT ${EXCHCXX_ENABLE_CUDA} )
    message( FATAL_ERROR "GauXC CUDA BINDINGS REQUIRE ExchCXX CUDA Bindings" )
  endif()

  if( ${GAUXC_ENABLE_SYCL} AND NOT ${EXCHCXX_ENABLE_SYCL} )
    message( FATAL_ERROR "GauXC SYCL BINDINGS REQUIRE ExchCXX SYCL Bindings" )
  endif()

endif()


