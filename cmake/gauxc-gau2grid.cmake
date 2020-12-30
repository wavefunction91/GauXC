if( GAUXC_ENABLE_GAU2GRID )

find_package( gau2grid CONFIG QUIET ) 
if( NOT gau2grid_FOUND )

  FetchContent_Declare(
    gau2grid
    GIT_REPOSITORY https://github.com/wavefunction91/gau2grid.git
    GIT_TAG        cmake-subproject
  )

  set( MAX_AM 6 CACHE STRING "" )
  set( DISABLE_PRAGMA ON CACHE BOOL "" )
  FetchContent_MakeAvailable( gau2grid )

endif()

if( NOT TARGET gau2grid::gg )
  message( STATUS "Something Went Horribly Wrong With Gau2Grid discovery!" )
endif()

endif() # If enabled
