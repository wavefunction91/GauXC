find_package( ExchCXX QUIET )
if( NOT ${ExchCXX_FOUND} )

  include( gauxc-dep-versions )

  # SYCL bindings only exist on the PR branch for now
  if( GAUXC_HAS_SYCL )
    set( GAUXC_EXCHCXX_REPOSITORY ${GAUXC_EXCHCXX_SYCL_REPOSITORY} )
    set( GAUXC_EXCHCXX_REVISION   ${GAUXC_EXCHCXX_SYCL_REVISION}   )
  endif()

  message( STATUS "Could not find ExchCXX... Building" )
  message( STATUS "EXCHCXX REPO = ${GAUXC_EXCHCXX_REPOSITORY}" )
  message( STATUS "EXCHCXX REV  = ${GAUXC_EXCHCXX_REVISION}"   )

  set( EXCHCXX_ENABLE_CUDA  ${GAUXC_HAS_CUDA} CACHE BOOL "" )
  set( EXCHCXX_ENABLE_HIP   ${GAUXC_HAS_HIP}  CACHE BOOL "" )
  set( EXCHCXX_ENABLE_SYCL  ${GAUXC_HAS_SYCL} CACHE BOOL "" )
  set( EXCHCXX_ENABLE_TESTS OFF               CACHE BOOL "" )

  if( GAUXC_HAS_SYCL AND GAUXC_SYCL_TARGET )
    set( EXCHCXX_SYCL_TARGET ${GAUXC_SYCL_TARGET} CACHE STRING "" )
  endif()

  FetchContent_Declare(
    exchcxx
    GIT_REPOSITORY ${GAUXC_EXCHCXX_REPOSITORY} 
    GIT_TAG        ${GAUXC_EXCHCXX_REVISION} 
  )

  FetchContent_MakeAvailable( exchcxx )


else()

  if( ${GAUXC_HAS_CUDA} AND NOT ${EXCHCXX_ENABLE_CUDA} )
    message( FATAL_ERROR "GauXC CUDA BINDINGS REQUIRE ExchCXX CUDA Bindings" )
  endif()

  if( ${GAUXC_HAS_HIP} AND NOT ${EXCHCXX_ENABLE_HIP} )
    message( FATAL_ERROR "GauXC HIP BINDINGS REQUIRE ExchCXX HIP Bindings" )
  endif()

  if( ${GAUXC_HAS_SYCL} AND NOT ${EXCHCXX_ENABLE_SYCL} )
    message( FATAL_ERROR "GauXC SYCL BINDINGS REQUIRE ExchCXX SYCL Bindings" )
  endif()

endif()


