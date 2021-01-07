cmake_minimum_required( VERSION 3.17 FATAL_ERROR )
project( gauxc_cmake_discovery LANGUAGES C CXX )

find_package( gauxc REQUIRED )
add_executable( gauxc_link_tester gauxc_link_tester.cxx )
target_link_libraries( gauxc_link_tester PUBLIC gauxc::gauxc )
