find_package( Eigen3 CONFIG QUIET )
if( NOT Eigen3_FOUND )
  
  FetchContent_Declare(
    eigen3
    URL https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz
  )

  FetchContent_GetPropertieS( eigen3 )
  if( NOT eigen3_POPULATED )
    FetchContent_Populate( eigen3 )
  endif()

  message( FATAL_ERROR "Eigen3 Pull Not Yet Configured" )

endif()

