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
include( FindPackageHandleStandardArgs )
include( CheckCXXCompilerFlag )
check_cxx_compiler_flag("-fsycl" CXX_HAS_FSYCL)

find_package_handle_standard_args( SYCL
  REQUIRED_VARS CXX_HAS_FSYCL
)

if( SYCL_FOUND AND NOT TARGET SYCL::SYCL )
  add_library( SYCL::SYCL INTERFACE IMPORTED )
  set_target_properties( SYCL::SYCL PROPERTIES
      INTERFACE_COMPILE_OPTIONS  "$<$<COMPILE_LANGUAGE:CXX>:-fsycl>"
      INTERFACE_LINK_OPTIONS     "-fsycl"
  )
endif()

# Flags every translation unit holding SYCL device code must be compiled with.
# GauXC applies these per source rather than target-wide; see gauxc-sycl.cmake.
set( SYCL_DEVICE_COMPILE_OPTIONS -fsycl )
