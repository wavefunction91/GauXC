#include <cmath>
#include <algorithm>

#include "cuda_runtime.h"

#include "cuda/cuda_device_properties.hpp"

namespace GauXC {
namespace cuda  {


uint32_t get_submat_cut_block(int32_t LDA, int32_t device) {
  int l2_cache_size;
  cudaDeviceGetAttribute(&l2_cache_size, cudaDevAttrL2CacheSize, device);

  int l2_block_size = (int) sqrt(0.75 * ((double) l2_cache_size / 8));
  int min_block_size = LDA / max_submat_blocks;

  int block_size = std::max(l2_block_size, min_block_size);
  block_size = std::min(block_size, LDA);

  return block_size;
}

uint32_t get_device_sm_count(int32_t device) {
  int num_sm;
  cudaDeviceGetAttribute(&num_sm, cudaDevAttrMultiProcessorCount, device);

  return num_sm;
}

}
}
