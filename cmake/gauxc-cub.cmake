if( GAUXC_HAS_CUDA )

  find_package( CUDAToolkit REQUIRED )
  if( CUDAToolkit_VERSION VERSION_LESS "11.0.0" )
    include( gauxc-dep-versions )

    message( STATUS "Building Local CUB Installation" )
    message( STATUS "CUB REPO = ${GAUXC_CUB_REPOSITORY}" )
    message( STATUS "CUB REV  = ${GAUXC_CUB_REVISION}"   )

    FetchContent_Declare(
      cub
      GIT_REPOSITORY ${GAUXC_CUB_REPOSITORY} 
      GIT_TAG        ${GAUXC_CUB_REVISION} 
    )

    FetchContent_GetProperties( cub )
    if( NOT cub_POPULATED )
      FetchContent_Populate( cub )
    endif()

    add_library( gauxc_cub INTERFACE IMPORTED )
    set_target_properties( gauxc_cub PROPERTIES 
      INTERFACE_INCLUDE_DIRECTORIES ${cub_SOURCE_DIR}
    )
  else()
    message( STATUS "Using CUB from CUDAToolkit" )
    message( STATUS "  CUDATOOLKIT VERSION = ${CUDAToolkit_VERSION}" )
  endif()

endif()
