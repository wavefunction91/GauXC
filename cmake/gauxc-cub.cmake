if( GAUXC_ENABLE_CUDA )

  FetchContent_Declare(
    cub
    GIT_REPOSITORY https://github.com/NVIDIA/cub.git
    GIT_TAG        1.10.0
  )

  FetchContent_GetProperties( cub )
  if( NOT cub_POPULATED )
    FetchContent_Populate( cub )
  endif()

  add_library( gauxc_cub INTERFACE IMPORTED )
  set_target_properties( gauxc_cub PROPERTIES 
    INTERFACE_INCLUDE_DIRECTORIES ${cub_SOURCE_DIR}
  )

endif()
