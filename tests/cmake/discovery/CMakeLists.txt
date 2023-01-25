#
# GauXC Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of
# any required approvals from the U.S. Dept. of Energy). All rights reserved.
#
# See LICENSE.txt for retails
#
cmake_minimum_required( VERSION 3.18 FATAL_ERROR )
project( gauxc_cmake_discovery LANGUAGES C CXX )

find_package( gauxc REQUIRED )
add_executable( gauxc_link_tester gauxc_link_tester.cxx )
target_link_libraries( gauxc_link_tester PUBLIC gauxc::gauxc )
