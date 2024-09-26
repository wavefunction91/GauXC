#
# GauXC Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of
# any required approvals from the U.S. Dept. of Energy). All rights reserved.
#
# See LICENSE.txt for details
#
set( GAUXC_OBARA_SAIKA_HOST_SRC
     src/integral_0.cxx
     src/integral_1.cxx
     src/integral_2.cxx
     src/integral_3.cxx
     src/integral_4.cxx
     src/integral_0_0.cxx
     src/integral_1_0.cxx
     src/integral_1_1.cxx
     src/integral_2_0.cxx
     src/integral_2_1.cxx
     src/integral_2_2.cxx
     src/integral_3_0.cxx
     src/integral_3_1.cxx
     src/integral_3_2.cxx
     src/integral_3_3.cxx
     src/integral_4_0.cxx
     src/integral_4_1.cxx
     src/integral_4_2.cxx
     src/integral_4_3.cxx
     src/integral_4_4.cxx
     src/obara_saika_integrals.cxx
     src/chebyshev_boys_computation.cxx
)

# This is required for compilation with NVHPC as it crashes with the O3 optimization flag
if (DEFINED GAUXC_OBARA_SAIKA_COMPILE_OPTIMIZATION_OPTIONS)
    foreach (flag "[\\/\\-]O3" "[\\/\\-]Ofast" "[\\/\\-]fast")
        string(REGEX REPLACE ${flag} ${GAUXC_OBARA_SAIKA_COMPILE_OPTIMIZATION_OPTIONS} CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
    endforeach()
    message("Setting Obara-Saika CMAKE_CXX_FLAGS_RELEASE to: ${CMAKE_CXX_FLAGS_RELEASE}")
endif()


target_sources( gauxc PRIVATE ${GAUXC_OBARA_SAIKA_HOST_SRC} )
target_include_directories( gauxc PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
)
