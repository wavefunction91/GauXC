if( GAUXC_ENABLE_GAU2GRID )
  if( NOT TARGET gau2grid::gg )
  
    # First try to find the package if target doesn't exist
    find_package( gau2grid CONFIG QUIET ) 
    
    if( NOT gau2grid_FOUND )
    
      message( STATUS "Could not find Gau2grid... Building" )
      
      if( GAUXC_FORCE_EXTERNAL_GAU2GRID )
        
        include( gauxc-dep-versions )
        
        message( STATUS "GAU2GRID URL = ${GAUXC_GAU2GRID_URL}" )
        
        FetchContent_Declare(
          gau2grid
          URL ${GAUXC_GAU2GRID_URL}
          URL_HASH SHA256=${GAUXC_GAU2GRID_SHA256}
          DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        )
        
        set( MAX_AM 6 CACHE STRING "" )
        set( DISABLE_PRAGMA ON CACHE BOOL "" )
        FetchContent_MakeAvailable( gau2grid )
        
        if( NOT TARGET gau2grid::gg )
          message( STATUS "Something Went Horribly Wrong With Gau2Grid discovery!" )
        endif()
      
      else()
      
        message( STATUS "Building Pregenerated Gau2grid" )
        add_subdirectory( ${PROJECT_SOURCE_DIR}/external/gau2grid ${PROJECT_BINARY_DIR}/external/gau2grid )
      
      endif()
    
    endif() # If not discoverable
  endif() # If target not present

  set(GAUXC_HAS_GAU2GRID TRUE CACHE BOOL "GauXC has Gau2Grid and will build bindings" FORCE)
endif() # If enabled
