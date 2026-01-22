# Find or download LibTorch
find_package(Torch QUIET)

if(NOT Torch_FOUND)
    message(STATUS "Torch not found. Downloading libtorch...")
   
    # Set libtorch version and download URL
    set(LIBTORCH_VERSION "2.9.1")
    set(USE_CUDA_LIBTORCH FALSE) #default is not to use the cuda version but cpu version
    if(GAUXC_HAS_CUDA)
        set(SUPPORTED_CUDA_VERSION_STRINGS "126" "128" "130")
        set(CUDA_VERSION_STRING ${CUDAToolkit_VERSION_MAJOR}${CUDAToolkit_VERSION_MINOR})
        if(CUDA_VERSION_STRING IN_LIST SUPPORTED_CUDA_VERSION_STRINGS)
	    set(USE_CUDA_LIBTORCH TRUE)
        else()
		message(WARNING "CUDA toolkit version is ${CUDAToolkit_VERSION}, for which there is no libtorch to download.")
		message(WARNING "Falling back to cpu version of libtorch.")
        endif()
    endif()
	
   
    # Determine the appropriate libtorch variant based on platform and CUDA availability
    if(UNIX AND NOT APPLE)
        # Linux
        if(USE_CUDA_LIBTORCH)
	    set(LIBTORCH_URL "https://download.pytorch.org/libtorch/cu${CUDA_VERSION_STRING}/libtorch-shared-with-deps-${LIBTORCH_VERSION}%2Bcu${CUDA_VERSION_STRING}.zip")
        else()
            set(LIBTORCH_URL "https://download.pytorch.org/libtorch/cpu/libtorch-shared-with-deps-${LIBTORCH_VERSION}%2Bcpu.zip")
        endif()
    elseif(APPLE)
        # macOS (CPU only)
        set(LIBTORCH_URL "https://download.pytorch.org/libtorch/cpu/libtorch-macos-arm64-${LIBTORCH_VERSION}.zip")
    elseif(WIN32)
        # Windows
        if(USE_CUDA_LIBTORCH)
	    set(LIBTORCH_URL "https://download.pytorch.org/libtorch/cu${CUDA_VERSION_STRING}/libtorch-win-shared-with-deps-${LIBTORCH_VERSION}%2Bcu${CUDA_VERSION_STRING}.zip")
        else()
            set(LIBTORCH_URL "https://download.pytorch.org/libtorch/cpu/libtorch-win-shared-with-deps-${LIBTORCH_VERSION}%2Bcpu.zip")
        endif()
    endif()
   
    set(LIBTORCH_DOWNLOAD_DIR "${CMAKE_BINARY_DIR}/libtorch")
    set(LIBTORCH_ZIP "${CMAKE_BINARY_DIR}/libtorch.zip")
   
    # Download libtorch
    if(NOT EXISTS ${LIBTORCH_DOWNLOAD_DIR})
        message(STATUS "Downloading libtorch from ${LIBTORCH_URL}")
        file(DOWNLOAD ${LIBTORCH_URL} ${LIBTORCH_ZIP}
            SHOW_PROGRESS
            STATUS DOWNLOAD_STATUS)
       
        list(GET DOWNLOAD_STATUS 0 STATUS_CODE)
        if(NOT STATUS_CODE EQUAL 0)
            message(FATAL_ERROR "Failed to download libtorch: ${DOWNLOAD_STATUS}")
        endif()
       
        # Extract libtorch
        message(STATUS "Extracting libtorch...")
        file(ARCHIVE_EXTRACT INPUT ${LIBTORCH_ZIP} DESTINATION ${CMAKE_BINARY_DIR})
       
        # Clean up zip file
        file(REMOVE ${LIBTORCH_ZIP})
    endif()
   
    # Set CMAKE_PREFIX_PATH to find the downloaded libtorch
    set(CMAKE_PREFIX_PATH "${LIBTORCH_DOWNLOAD_DIR};${CMAKE_PREFIX_PATH}")
   
    # Find Torch package again
    find_package(Torch REQUIRED)
endif()

message(STATUS "Using Torch version: ${Torch_VERSION}")
