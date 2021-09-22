#pragma once
#include "device/device_queue.hpp"

namespace GauXC {

void symmetrize_matrix( int32_t N, double* A, size_t LDA, device_queue queue );


}
