#pragma once
#include "hip/hip_runtime.h"

namespace GauXC {

void symmetrize_matrix( int32_t N, double* A, size_t LDA, hipStream_t stream );


}
