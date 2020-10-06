if( GAUXC_ENABLE_GAU2GRID )

add_library( gauxc_gau2grid INTERFACE )

find_package( gau2grid CONFIG QUIET ) 
if( NOT gau2grid_FOUND )

  FetchContent_Declare(
    gau2grid
    GIT_REPOSITORY https://github.com/dgasmith/gau2grid.git
  )

  set( MAX_AM 6 CACHE STRING "" )
  FetchContent_GetProperties( gau2grid )
  if( NOT gau2grid_POPULATED )
    FetchContent_Populate( gau2grid )
    file( READ ${gau2grid_SOURCE_DIR}/CMakeLists.txt GAU2GRID_CMAKE_LISTS )
    string( REPLACE "CMAKE_SOURCE_DIR" "PROJECT_SOURCE_DIR" GAU2GRID_CORRECT "${GAU2GRID_CMAKE_LISTS}" )
    file( WRITE ${gau2grid_SOURCE_DIR}/CMakeLists.txt "${GAU2GRID_CORRECT}" )
    add_subdirectory( ${gau2grid_SOURCE_DIR} ${gau2grid_BINARY_DIR} )
    target_compile_definitions( gg PRIVATE $<BUILD_INTERFACE:__GG_NO_PRAGMA> )
  endif()



  target_link_libraries( gauxc_gau2grid INTERFACE gg )
  target_include_directories( gauxc_gau2grid INTERFACE $<BUILD_INTERFACE:${gau2grid_BINARY_DIR}> )

else()

  target_link_libraries( gauxc_gau2grid INTERFACE gau2grid::gg )

endif()

endif() # If enabled
