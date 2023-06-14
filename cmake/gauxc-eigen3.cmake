find_package( Eigen3 CONFIG HINTS ${EIGEN3_ROOT_DIR} )
if( NOT Eigen3_FOUND )
  
  message( STATUS "Could Not Find Eigen3... Building" )
  message( STATUS "EIGEN3 REPO = https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz" )
  #message( STATUS "EIGEN3 REV  = "   )

  FetchContent_Declare(
    eigen3
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
  )

  FetchContent_GetProperties( eigen3 )
  if( NOT eigen3_POPULATED )
    FetchContent_Populate( eigen3 )
  endif()

  #message( FATAL_ERROR "Eigen3 Pull Not Yet Configured" )
  add_library( Eigen3::Eigen INTERFACE IMPORTED )
  set_target_properties( Eigen3::Eigen PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${eigen3_SOURCE_DIR}
  )

endif()

