find_package( IntegratorXX QUIET )
if( NOT ${IntegratorXX_FOUND} )

  include( gauxc-dep-versions )

  message( STATUS "Could not find IntegratorXX... Building" )
  message( STATUS "INTEGRATORXX REPO = ${GAUXC_INTEGRATORXX_REPOSITORY}" )
  message( STATUS "INTEGRATORXX REV  = ${GAUXC_INTEGRATORXX_REVISION}"   )

  set( INTEGRATORXX_ENABLE_TESTS OFF CACHE BOOL "" )
  FetchContent_Declare(
    integratorxx
    GIT_REPOSITORY ${GAUXC_INTEGRATORXX_REPOSITORY} 
    GIT_TAG        ${GAUXC_INTEGRATORXX_REVISION} 
  )

  FetchContent_MakeAvailable( integratorxx )

endif()


