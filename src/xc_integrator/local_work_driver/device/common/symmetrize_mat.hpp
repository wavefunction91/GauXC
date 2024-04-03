/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include "device/device_queue.hpp"

namespace GauXC {

void symmetrize_matrix( int32_t N, double* A, size_t LDA, device_queue queue );
void symmetrize_matrix_inc( int32_t N, double* A, size_t LDA, device_queue queue );


}
