include( FetchContent )
include( gauxc-dep-versions )
FetchContent_Declare( linalg-cmake-modules 
  GIT_REPOSITORY ${GAUXC_LINALG_MODULES_REPOSITORY} 
  GIT_TAG        ${GAUXC_LINALG_MODULES_REVISION} 
)
FetchContent_GetProperties( linalg-cmake-modules )
if( NOT linalg-cmake-modules_POPULATED )
  FetchContent_Populate( linalg-cmake-modules )
  list( PREPEND CMAKE_MODULE_PATH ${linalg-cmake-modules_SOURCE_DIR} )
endif()
