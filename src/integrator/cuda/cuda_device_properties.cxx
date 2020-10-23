#include <cmath>
#include <algorithm>

#include "cuda_runtime.h"

#include "cuda_device_properties.hpp"

namespace GauXC {
namespace cuda  {


uint32_t get_submat_cut_block(int32_t LDA) {
  int l2_cache_size;
  cudaDeviceGetAttribute(&l2_cache_size, cudaDevAttrL2CacheSize, 0);

  int l2_block_size = (int) sqrt(0.75 * ((double) l2_cache_size / 8));
  int min_block_size = LDA / max_submat_blocks;

  int block_size = std::max(l2_block_size, min_block_size);
  block_size = std::min(block_size, LDA);

  return block_size;
}

}
}
