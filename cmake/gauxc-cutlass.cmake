# Check that only CUDA CC 8.0+ is enabled
foreach( cuda_arch ${CMAKE_CUDA_ARCHITECTURES} )
  if( NOT cuda_arch GREATER_EQUAL 80 )
    message(FATAL_ERROR "GauXC Requires CUDA CC >= 8.0 For CUTLASS")
  endif()
endforeach()

include( gauxc-dep-versions )

message( STATUS "Building Local CUTLASS Installation" )
message( STATUS "CUTLASS REPO = ${GAUXC_CUTLASS_REPOSITORY}" )
message( STATUS "CUTLASS REV  = ${GAUXC_CUTLASS_REVISION}"   )

FetchContent_Declare(
  cutlass
  GIT_REPOSITORY ${GAUXC_CUTLASS_REPOSITORY} 
  GIT_TAG        ${GAUXC_CUTLASS_REVISION} 
)

FetchContent_GetProperties( cutlass )
if( NOT cutlass_POPULATED )
  FetchContent_Populate( cutlass )
endif()



add_library( gauxc_cutlass INTERFACE IMPORTED )
set_target_properties( gauxc_cutlass PROPERTIES 
  INTERFACE_INCLUDE_DIRECTORIES 
    "${cutlass_SOURCE_DIR}/include;${cutlass_SOURCE_DIR}/tools/util/include"
)

set(GAUXC_HAS_CUTLASS TRUE CACHE BOOL "GauXC has CUTLASS" FORCE) 
