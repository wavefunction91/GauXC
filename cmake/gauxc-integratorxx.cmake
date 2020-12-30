find_package( IntegratorXX QUIET )
if( NOT ${IntegratorXX_FOUND} )

  set( INTEGRATORXX_ENABLE_TESTS OFF CACHE BOOL "" )
  FetchContent_Declare(
    integratorxx
    GIT_REPOSITORY 
      https://github.com/wavefunction91/IntegratorXX.git
    GIT_TAG 
      cmake-fixup
  )

  FetchContent_MakeAvailable( integratorxx )

endif()


