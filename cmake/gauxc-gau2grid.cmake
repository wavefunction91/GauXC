if( GAUXC_ENABLE_GAU2GRID )

if( GAUXC_FORCE_EXTERNAL_GAU2GRID )

  find_package( gau2grid CONFIG QUIET ) 
  if( NOT gau2grid_FOUND )
  
    include( gauxc-dep-versions )
  
    message( STATUS "Could not find Gau2grid... Building" )
    message( STATUS "GAU2GRID REPO = ${GAUXC_GAU2GRID_REPOSITORY}" )
    message( STATUS "GAU2GRID REV  = ${GAUXC_GAU2GRID_REVISION}"   )
  
    FetchContent_Declare(
      gau2grid
      GIT_REPOSITORY ${GAUXC_GAU2GRID_REPOSITORY} 
      GIT_TAG        ${GAUXC_GAU2GRID_REVISION} 
    )
  
    set( MAX_AM 6 CACHE STRING "" )
    set( DISABLE_PRAGMA ON CACHE BOOL "" )
    FetchContent_MakeAvailable( gau2grid )
  
  endif()
  
  if( NOT TARGET gau2grid::gg )
    message( STATUS "Something Went Horribly Wrong With Gau2Grid discovery!" )
  endif()

else()

  add_subdirectory( ${PROJECT_SOURCE_DIR}/external/gau2grid ${PROJECT_BINARY_DIR}/external/gau2grid )

endif()

endif() # If enabled
