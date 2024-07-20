#
# GauXC Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of
# any required approvals from the U.S. Dept. of Energy). All rights reserved.
#
# See LICENSE.txt for details
#

list (PREPEND CMAKE_PREFIX_PATH ${ROCM_PATH} ${ROCM_PATH}/hip ${ROCM_PATH}/hipblas)
set(GPU_TARGETS "${CMAKE_HIP_ARCHITECTURES}" CACHE STRING "AMD GPU targets to compile for")

find_package( hip REQUIRED )
find_package( hipblas REQUIRED )

list(REMOVE_AT CMAKE_PREFIX_PATH 0 1 2)

target_sources( gauxc PRIVATE hip_backend.cxx )
target_link_libraries( gauxc PUBLIC hip::host roc::hipblas )
