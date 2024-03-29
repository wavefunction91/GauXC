find_package( cereal QUIET )
if( NOT cereal_FOUND )

  include( gauxc-dep-versions )

  message( STATUS "Could not find Cereal... Building" )
  message( STATUS "CEREAL REPO = ${GAUXC_CEREAL_REPOSITORY}" )
  message( STATUS "CEREAL REV  = ${GAUXC_CEREAL_REVISION}"   )

  FetchContent_Declare(
    cereal
    GIT_REPOSITORY ${GAUXC_CEREAL_REPOSITORY} 
    GIT_TAG        ${GAUXC_CEREAL_REVISION} 
  )

  FetchContent_GetProperties(cereal)
  if(NOT cereal_POPULATED)
    FetchContent_Populate( cereal )
    add_library( cereal INTERFACE IMPORTED )
    set_target_properties( cereal PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${cereal_SOURCE_DIR}/include"
      INTERFACE_COMPILE_DEFINITIONS "CEREAL_THREAD_SAFE=1;GAUXC_HAS_CEREAL=1"
    )
  endif()

else()

  target_compile_definitions( cereal INTERFACE
    "CEREAL_THREAD_SAFE=1;GAUXC_HAS_CEREAL=1" 
  )

endif()
