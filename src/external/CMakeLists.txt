#
# GauXC Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of
# any required approvals from the U.S. Dept. of Energy).
#
# (c) 2024-2025, Microsoft Corporation
#
# All rights reserved.
#
# See LICENSE.txt for details
#
if( GAUXC_ENABLE_HDF5 )
  include(FetchContent)
  find_package(HDF5)
  if(HDF5_FOUND)
    set(GAUXC_HAS_HDF5 TRUE CACHE BOOL "" FORCE)
    message(STATUS "Enabling HDF5 Bindings")
    find_package(HighFive QUIET)
    if(NOT HighFive_FOUND)
      message(STATUS "HighFive REPO = ${GAUXC_HIGHFIVE_REPOSITORY}")
      message(STATUS "HighFive REV  = ${GAUXC_HIGHFIVE_REVISION}  ")
      FetchContent_Declare( HighFive
        GIT_REPOSITORY ${GAUXC_HIGHFIVE_REPOSITORY}
        GIT_TAG        ${GAUXC_HIGHFIVE_REVISION}
      )

      set(HIGHFIVE_USE_BOOST OFF CACHE BOOL "" )
      set(HIGHFIVE_UNIT_TESTS OFF CACHE BOOL "" )
      set(HIGHFIVE_EXAMPLES OFF CACHE BOOL "" )
      #set(HIGHFIVE_PARALLEL_HDF5 ON CACHE BOOL "" )
      set(HIGHFIVE_BUILD_DOCS OFF CACHE BOOL "" )
      FetchContent_MakeAvailable( HighFive )

      # HighFive propagates HDF5_DEFINITIONS via its libdeps target.
      # CMake's FindHDF5 module sets H5_BUILT_AS_DYNAMIC_LIB on Windows
      # when HDF5_USE_STATIC_LIBRARIES is not set, even for static libs.
      # Correct this when linking statically.
      if(WIN32 AND HDF5_PROVIDES_STATIC_LIBS AND TARGET libdeps)
        get_target_property(_libdeps_defs libdeps INTERFACE_COMPILE_DEFINITIONS)
        if(_libdeps_defs)
          list(REMOVE_ITEM _libdeps_defs "H5_BUILT_AS_DYNAMIC_LIB")
          list(APPEND _libdeps_defs "H5_BUILT_AS_STATIC_LIB")
          set_target_properties(libdeps PROPERTIES INTERFACE_COMPILE_DEFINITIONS "${_libdeps_defs}")
        endif()
      endif()

    endif()
    target_sources( gauxc PRIVATE hdf5_write.cxx hdf5_read.cxx )
    if( GAUXC_ENABLE_C )
      target_sources( gauxc PRIVATE c_hdf5_read.cxx c_hdf5_write.cxx )
    endif()
    if( GAUXC_ENABLE_FORTRAN )
      target_sources( gauxc PRIVATE hdf5_read.f90 hdf5_write.f90 )
    endif()
    target_link_libraries( gauxc PUBLIC HighFive )

    # When linking HDF5 statically on Windows, HDF5's transitive
    # dependencies (zlib, szip/aec, shlwapi) must be linked explicitly.
    if(WIN32 AND HDF5_PROVIDES_STATIC_LIBS)
      find_library(ZLIB_LIBRARY NAMES zlib z)
      find_library(SZIP_LIBRARY NAMES szip-static szip sz)
      find_library(AEC_LIBRARY  NAMES aec-static aec)
      target_link_libraries( gauxc PUBLIC
        ${ZLIB_LIBRARY} ${SZIP_LIBRARY} ${AEC_LIBRARY} shlwapi )
    endif()
  else()
    message(WARNING "GAUXC_ENABLE_HDF5 was enabled, but HDF5 was not found, Disabling HDF5 Bindings")
  endif()
else()
  message(STATUS "Disabling HDF5 Bindings")
endif()
