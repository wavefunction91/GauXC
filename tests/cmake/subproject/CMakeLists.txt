cmake_minimum_required( VERSION 3.17 FATAL_ERROR )
project( gauxc_cmake_subproject LANGUAGES C CXX )


if( NOT GIT_REVISION )
  # Get current Git Revision
  find_package( Git REQUIRED )
  execute_process( 
    COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
    OUTPUT_VARIABLE GIT_REVISION
    ERROR_QUIET
  )
  string( STRIP "${GIT_REVISION}" GIT_REVISION )
endif()

message( STATUS "Pulling GIT_REVISION = ${GIT_REVISION}" )

include( FetchContent )
FetchContent_Declare( gauxc
  GIT_REPOSITORY https://github.com/wavefunction91/gauxc.git
  GIT_TAG        ${GIT_REVISION}
)
FetchContent_MakeAvailable( gauxc )

add_executable( gauxc_link_tester gauxc_link_tester.cxx )
target_link_libraries( gauxc_link_tester PUBLIC gauxc::gauxc )
