#
# GauXC Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of
# any required approvals from the U.S. Dept. of Energy). All rights reserved.
#
# See LICENSE.txt for retails
#
if( GAUXC_ENABLE_HDF5 )
  include(FetchContent)
  find_package(HDF5 REQUIRED)
  FetchContent_Declare( HighFive
    GIT_REPOSITORY https://github.com/BlueBrain/HighFive.git
    GIT_TAG        805f0e13d09b47c4b01d40682621904aa3b31bb8
  )
  
  set(HIGHFIVE_USE_BOOST OFF CACHE BOOL "" )
  set(HIGHFIVE_UNIT_TESTS OFF CACHE BOOL "" )
  set(HIGHFIVE_EXAMPLES OFF CACHE BOOL "" )
  #set(HIGHFIVE_PARALLEL_HDF5 ON CACHE BOOL "" )
  set(HIGHFIVE_BUILD_DOCS OFF CACHE BOOL "" )
  FetchContent_MakeAvailable( HighFive )
  
  target_sources( gauxc PRIVATE hdf5_write.cxx hdf5_read.cxx )
  target_link_libraries( gauxc PUBLIC HighFive )
endif()
