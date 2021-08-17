target_sources( gauxc PRIVATE 
  replicated_xc_integrator_impl.cxx 
)

add_subdirectory(host)

if(GAUXC_ENABLE_DEVICE)
  add_subdirectory(device)
endif()

target_include_directories( gauxc
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
)
