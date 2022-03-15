# This CMake harness is meant for use with the GauXC library
# and is released under the terms of the 3-clause BSD license

target_sources( gauxc PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/generated_source/gau2grid_phi.c
  ${CMAKE_CURRENT_SOURCE_DIR}/generated_source/gau2grid_orbital.c
  ${CMAKE_CURRENT_SOURCE_DIR}/generated_source/gau2grid_deriv1.c
  ${CMAKE_CURRENT_SOURCE_DIR}/generated_source/gau2grid_deriv2.c
  ${CMAKE_CURRENT_SOURCE_DIR}/generated_source/gau2grid_deriv3.c
  ${CMAKE_CURRENT_SOURCE_DIR}/generated_source/gau2grid_transform.c
  ${CMAKE_CURRENT_SOURCE_DIR}/generated_source/gau2grid_helper.c )


target_compile_definitions( gauxc PRIVATE $<BUILD_INTERFACE:__GG_NO_PRAGMA> )
target_include_directories( gauxc 
  PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/generated_source>
)
