# Find or download LibTorch
find_package(Torch QUIET)

if(NOT Torch_FOUND)
    message(STATUS "Torch not found. Downloading libtorch...")
   
    # Set libtorch version and download URL
    set(LIBTORCH_VERSION "2.9.1")
   
    # Determine the appropriate libtorch variant based on platform and CUDA availability
    if(UNIX AND NOT APPLE)
        # Linux
        if(GAUXC_HAS_CUDA)
            set(LIBTORCH_URL "https://download.pytorch.org/libtorch/cu121/libtorch-cxx11-abi-shared-with-deps-${LIBTORCH_VERSION}%2Bcu121.zip")
        else()
            set(LIBTORCH_URL "https://download.pytorch.org/libtorch/cpu/libtorch-shared-with-deps-${LIBTORCH_VERSION}%2Bcpu.zip"
        endif()
    elseif(APPLE)
        # macOS (CPU only)
        set(LIBTORCH_URL "https://download.pytorch.org/libtorch/cpu/libtorch-macos-arm64-${LIBTORCH_VERSION}.zip")
    elseif(WIN32)
        # Windows
        if(GAUXC_HAS_CUDA)
            set(LIBTORCH_URL "https://download.pytorch.org/libtorch/cu121/libtorch-win-shared-with-deps-${LIBTORCH_VERSION}%2Bcu121.zip")
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
