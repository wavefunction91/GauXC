# Check that only CUDA CC 8.0+ is enabled
foreach( cuda_arch ${CMAKE_CUDA_ARCHITECTURES} )
  if( NOT cuda_arch GREATER_EQUAL 80 )
    message(FATAL_ERROR "GauXC Requires CUDA CC >= 8.0 For CUTLASS")
  endif()
endforeach()

include( gauxc-dep-versions )

message( STATUS "Building Local CUTLASS Installation" )
message( STATUS "CUTLASS URL = ${GAUXC_CUTLASS_URL}" )

FetchContent_Declare(
  cutlass
  URL ${GAUXC_CUTLASS_URL}
  URL_HASH SHA256=${GAUXC_CUTLASS_SHA256}
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
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
