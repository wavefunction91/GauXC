cmake_minimum_required( VERSION 3.17 FATAL_ERROR )

include(FetchContent)
project( GauXC VERSION 0.0.1 LANGUAGES C CXX )

option( GAUXC_ENABLE_CUDA "Enable CUDA Bindings" OFF )
option( GAUXC_ENABLE_TESTS "Enable Unit Tests"   ON  )

if( GAUXC_ENABLE_CUDA )
  enable_language( CUDA )
endif()

add_subdirectory( src )

if( GAUXC_ENABLE_TESTS )
  enable_testing()
  add_subdirectory( tests )
endif()
