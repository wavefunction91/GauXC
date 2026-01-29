find_package( Eigen3 CONFIG HINTS ${EIGEN3_ROOT_DIR} )
if( NOT Eigen3_FOUND )
  
  message( STATUS "Could Not Find Eigen3... Building" )
  message( STATUS "EIGEN3 URL = ${GAUXC_EIGEN3_URL}" )
  #message( STATUS "EIGEN3 REV  = "   )

  FetchContent_Declare(
    eigen3
    URL ${GAUXC_EIGEN3_URL}
    URL_HASH SHA256=${GAUXC_EIGEN3_SHA256}
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
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

