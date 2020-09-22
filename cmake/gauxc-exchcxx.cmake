find_package( ExchCXX QUIET )
if( NOT ${ExchCXX_FOUND} )

  set( EXCHCXX_ENABLE_CUDA  ${GAUXC_ENABLE_CUDA} CACHE BOOL "" )
  set( EXCHCXX_ENABLE_SYCL  ${GAUXC_ENABLE_SYCL} CACHE BOOL "" )
  set( EXCHCXX_ENABLE_HIP   ${GAUXC_ENABLE_HIP}  CACHE BOOL "" )
  set( EXCHCXX_ENABLE_TESTS OFF                  CACHE BOOL "" )

  FetchContent_Declare(
    exchcxx
    GIT_REPOSITORY 
      https://github.com/wavefunction91/ExchCXX.git
    GIT_TAG 
      perf_port
  )

  FetchContent_MakeAvailable( exchcxx )


else()

  if( ${GAUXC_ENABLE_CUDA} AND NOT ${EXCHCXX_ENABLE_CUDA} )
    message( FATAL_ERROR "GauXC CUDA BINDINGS REQUIRE ExchCXX CUDA Bindings" )
  endif()

  if( ${GAUXC_ENABLE_SYCL} AND NOT ${EXCHCXX_ENABLE_SYCL} )
    message( FATAL_ERROR "GauXC SYCL BINDINGS REQUIRE ExchCXX SYCL Bindings" )
  endif()

  if( ${GAUXC_ENABLE_HIP} AND NOT ${EXCHCXX_ENABLE_HIP} )
    message( FATAL_ERROR "GauXC HIP BINDINGS REQUIRE ExchCXX HIP Bindings" )
  endif()

endif()


