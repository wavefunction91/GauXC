find_package( cereal QUIET )
if( NOT cereal_FOUND )

  FetchContent_Declare(
    cereal
    GIT_REPOSITORY https://github.com/USCiLab/cereal.git
  )
  FetchContent_GetProperties(cereal)
  if(NOT cereal_POPULATED)
    FetchContent_Populate( cereal )
    add_library( cereal INTERFACE IMPORTED )
    set_target_properties( cereal PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${cereal_SOURCE_DIR}/include"
      INTERFACE_COMPILE_DEFINITIONS "CEREAL_THREAD_SAFE=1;GAUXC_ENABLE_CEREAL=1"
    )
  endif()

else()

  target_compile_definitions( cereal INTERFACE
    "CEREAL_THREAD_SAFE=1;GAUXC_ENABLE_CEREAL=1" 
  )

endif()
