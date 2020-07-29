find_package( ExchCXX QUIET )
if( NOT ${ExchCXX_FOUND} )

  set( EXCHCXX_ENABLE_CUDA  ${GAUXC_ENABLE_CUDA} CACHE BOOL "" )
  set( EXCHCXX_ENABLE_TESTS OFF                  CACHE BOOL "" )

  FetchContent_Declare(
    exchcxx
    GIT_REPOSITORY 
      https://github.com/wavefunction91/ExchCXX.git
    GIT_TAG 
      master
  )

  FetchContent_MakeAvailable( exchcxx )


else()

  if( ${GAUXC_ENABLE_CUDA} AND NOT ${EXCHCXX_ENABLE_CUDA} )
    message( FATAL_ERROR "GauXC CUDA BINDINGS REQUIRE ExchCXX CUDA Bindings" )
  endif()

endif()


