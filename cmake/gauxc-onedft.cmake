if(NOT TARGET nlohmann_json::nlohmann_json)
  find_package(nlohmann_json)
  if( NOT nlohmann_json_FOUND )

    message( STATUS "Could Not Find nlohmann_json... Building" )
    message( STATUS "NLOHMANN_JSON URL = ${GAUXC_NLOHMANN_JSON_URL}" )

    FetchContent_Declare(
      nlohmann_json
      URL ${GAUXC_NLOHMANN_JSON_URL}
      URL_HASH SHA256=${GAUXC_NLOHMANN_JSON_SHA256}
      DOWNLOAD_EXTRACT_TIMESTAMP ON
    )

    FetchContent_GetProperties( nlohmann_json )
    if( NOT nlohmann_json_POPULATED )
      FetchContent_Populate( nlohmann_json )
    endif()

    add_library( nlohmann_json::nlohmann_json INTERFACE IMPORTED )
    set_target_properties( nlohmann_json::nlohmann_json PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES ${nlohmann_json_SOURCE_DIR}/include
    )
  endif()
endif()

# store and restore CMAKE_CUDA_ARCHITECTURES and CMAKE_CUDA_FLAGS if Torch clobbers them
set(_PREV_CUDA_ARCHS "${CMAKE_CUDA_ARCHITECTURES}")
set(_PREV_CUDA_FLAGS "${CMAKE_CUDA_FLAGS}")
find_package(Torch REQUIRED)
# Restore CMAKE_CUDA_ARCHITECTURES (Torch may set it to OFF)
if(NOT "${CMAKE_CUDA_ARCHITECTURES}" STREQUAL "${_PREV_CUDA_ARCHS}")
  set(CMAKE_CUDA_ARCHITECTURES "${_PREV_CUDA_ARCHS}" CACHE STRING "" FORCE)
  message(WARNING "Torch changed CMAKE_CUDA_ARCHITECTURES. Restored previous value: ${CMAKE_CUDA_ARCHITECTURES}")
endif()
# Strip Torch-injected -gencode flags from CMAKE_CUDA_FLAGS (PyTorch issue #71379)
string(REGEX REPLACE " -gencode [^ ]+" "" _cleaned_cuda_flags "${CMAKE_CUDA_FLAGS}")
if(NOT "${_cleaned_cuda_flags}" STREQUAL "${CMAKE_CUDA_FLAGS}")
  set(CMAKE_CUDA_FLAGS "${_cleaned_cuda_flags}" CACHE STRING "" FORCE)
  message(WARNING "Stripped Torch-injected -gencode flags from CMAKE_CUDA_FLAGS")
endif()
list(REMOVE_ITEM TORCH_LIBRARIES torch::nvtoolsext)
message(STATUS "Torch libraries without nvtoolsext: ${TORCH_LIBRARIES}")
