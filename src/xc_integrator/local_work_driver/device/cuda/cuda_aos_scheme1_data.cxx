#include "cuda_aos_scheme1.hpp"
#include "buffer_adaptor.hpp"
#include <gauxc/util/div_ceil.hpp>
#include "device/cuda/cuda_backend.hpp"

namespace GauXC {

CudaAoSScheme1::Data::~Data() noexcept = default;

CudaAoSScheme1::Data::Data(bool batch_l3_blas) :
  Scheme1DataBase( std::make_unique<CUDABackend>(), batch_l3_blas ) { }

size_t CudaAoSScheme1::Data::get_ldatoms() {
  constexpr auto weight_unroll = CudaAoSScheme1::weight_unroll;
  return util::div_ceil( global_dims.natoms, weight_unroll ) * weight_unroll;
}

size_t CudaAoSScheme1::Data::get_rab_align() {
  return sizeof(double2);
}


size_t CudaAoSScheme1::Data::get_submat_chunk_size(int32_t LDA, int32_t dev_id) {

  constexpr auto max_submat_blocks = CudaAoSScheme1::max_submat_blocks;

  int l2_cache_size;
  cudaDeviceGetAttribute(&l2_cache_size, cudaDevAttrL2CacheSize, dev_id);

  int l2_block_size = (int) sqrt(0.75 * ((double) l2_cache_size / 8));
  int min_block_size = LDA / max_submat_blocks;

  int block_size = std::max(l2_block_size, min_block_size);
  block_size = std::min(block_size, LDA);

  return block_size;

}


}
