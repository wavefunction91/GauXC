find_package(nlohmann_json REQUIRED)
# store and restore CMAKE_CUDA_ARCHITECTURES if Torch clobbers it

set(_PREV_CUDA_ARCHS "${CMAKE_CUDA_ARCHITECTURES}")
find_package(Torch REQUIRED)
if(CMAKE_CUDA_ARCHITECTURES STREQUAL "OFF")
  set(CMAKE_CUDA_ARCHITECTURES "${_PREV_CUDA_ARCHS}" CACHE STRING "Restore CUDA archs after Torch override" FORCE)
  message(WARNING "Torch set CMAKE_CUDA_ARCHITECTURES to OFF. Restored previous value: ${CMAKE_CUDA_ARCHITECTURES}")
endif()
list(REMOVE_ITEM TORCH_LIBRARIES torch::nvtoolsext)
message(STATUS "Torch libraries without nvtoolsext: ${TORCH_LIBRARIES}")
