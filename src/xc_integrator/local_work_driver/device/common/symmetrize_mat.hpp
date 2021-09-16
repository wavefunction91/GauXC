#pragma once
#include "device/type_erased_queue.hpp"

namespace GauXC {

void symmetrize_matrix( int32_t N, double* A, size_t LDA, type_erased_queue queue );


}
