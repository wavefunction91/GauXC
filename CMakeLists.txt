cmake_minimum_required( VERSION 3.18 FATAL_ERROR )

include(FetchContent)
set( FETCHCONTENT_UPDATES_DISCONNECTED ON CACHE BOOL "Disable FC Updates" )

project( GauXC VERSION 0.0.1 LANGUAGES C CXX )

# Place local modules in the path
list( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake )
list( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules )
include( gauxc-linalg-modules )

# Populate BUILD_TESTING prior to dependencies to avoid clash
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  include(CTest)
endif()



# GauXC Options
option( GAUXC_ENABLE_HOST     "Enable Host Integrator"      ON  )
option( GAUXC_ENABLE_CUDA     "Enable CUDA Bindings"        OFF )
option( GAUXC_ENABLE_MAGMA    "Enable MAGMA Linear Algebra" ON  )
option( GAUXC_ENABLE_NCCL     "Enable NCCL Collectives    " OFF  )
option( GAUXC_ENABLE_MPI      "Enable MPI Bindings"         ON  )
option( GAUXC_ENABLE_TESTS    "Enable Unit Tests"           ON  )
option( GAUXC_ENABLE_GAU2GRID "Enable Gau2Grid Collocation" ON  )

if( NOT GAUXC_ENABLE_GAU2GRID )
  message( FATAL_ERROR "Gau2Grid is currently a required dependency which 
                        will be made optional in a future releast of GauXC [WIP]" ) 
endif()

# Decided if we're compiling device bindings
if( GAUXC_ENABLE_CUDA )
  set( GAUXC_ENABLE_DEVICE TRUE CACHE BOOL "Enable Device Code" )
else()
  set( GAUXC_ENABLE_DEVICE FALSE CACHE BOOL "Enable Device Code" )
endif()



if( NOT (${GAUXC_ENABLE_HOST} OR ${GAUXC_ENABLE_DEVICE}) )
  message( FATAL_ERROR "Neither Host nor Device Integrators have been enabled!" )
endif()

if( GAUXC_ENABLE_CUDA )
  enable_language( CUDA )
endif()

add_subdirectory( src )

if( CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND GAUXC_ENABLE_TESTS AND BUILD_TESTING )
  add_subdirectory( tests )
endif()
