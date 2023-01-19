if(GAUXC_ENABLE_DEVICE)
  add_subdirectory( device )
  target_include_directories( gauxc
    PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
endif()
