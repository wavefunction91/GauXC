#pragma once

namespace GauXC {

void symmetrize_matrix( int32_t N, double* A, size_t LDA, cudaStream_t stream );


}
