add_library( gauxc_catch2 INTERFACE )

find_package( Catch2 CONFIG QUIET )
if( NOT Catch2_FOUND )

  FetchContent_Declare(
    catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  )
  
  set(CATCH_BUILD_TESTING OFF CACHE BOOL "Build SelfTest project" FORCE)
  set(CATCH_INSTALL_DOCS OFF CACHE BOOL "Install documentation alongside library" FORCE)
  set(CATCH_INSTALL_HELPERS OFF CACHE BOOL "Install contrib alongside library" FORCE)

  FetchContent_MakeAvailable( catch2 )
  target_link_libraries( gauxc_catch2 INTERFACE Catch2::Catch2 )

else()

  target_link_libraries( gauxc_catch2 INTERFACE Catch2::Catch2 )

endif()

#add_library( catch2_main STATIC ut_main.cxx)
#target_link_libraries( catch2_main PUBLIC gauxc_catch2 )
#
#add_executable( gauxc_moltypes_test moltypes_test.cxx )
#target_link_libraries( gauxc_moltypes_test PUBLIC gauxc catch2_main )
#
#add_executable( gauxc_grid_test grid_test.cxx )
#target_link_libraries( gauxc_grid_test PUBLIC gauxc catch2_main )

add_executable( gauxc_test ut_main.cxx moltypes_test.cxx grid_test.cxx molgrid_test.cxx basisset_test.cxx load_balancer_test.cxx xc_integrator.cxx standards.cxx )
target_link_libraries( gauxc_test PUBLIC gauxc gauxc_catch2 )


set( GAUXC_REF_DATA_PATH "${PROJECT_SOURCE_DIR}/tests/ref_data" )
configure_file( ut_common.hpp.in ${PROJECT_BINARY_DIR}/tests/ut_common.hpp )
target_include_directories( gauxc_test PRIVATE ${PROJECT_BINARY_DIR}/tests )
target_include_directories( gauxc_test PRIVATE ${PROJECT_SOURCE_DIR}/tests )
