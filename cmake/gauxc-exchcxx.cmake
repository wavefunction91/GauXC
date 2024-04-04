find_package( ExchCXX QUIET )
if( NOT ${ExchCXX_FOUND} )

  include( gauxc-dep-versions )

  message( STATUS "Could not find ExchCXX... Building" )
  message( STATUS "EXCHCXX REPO = ${GAUXC_EXCHCXX_REPOSITORY}" )
  message( STATUS "EXCHCXX REV  = ${GAUXC_EXCHCXX_REVISION}"   )

  set( EXCHCXX_ENABLE_CUDA  ${GAUXC_HAS_CUDA} CACHE BOOL "" )
  set( EXCHCXX_ENABLE_HIP   ${GAUXC_HAS_HIP}  CACHE BOOL "" )
  set( EXCHCXX_ENABLE_TESTS OFF               CACHE BOOL "" )

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

endif()


