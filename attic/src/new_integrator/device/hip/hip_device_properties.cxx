#include <cmath>
#include <algorithm>

#include "hip_runtime.h"

#include "device/hip/hip_device_properties.hpp"

namespace GauXC {
namespace hip  {


uint32_t get_submat_cut_block(int32_t LDA, int32_t device) {
  int l2_cache_size;
  hipDeviceGetAttribute(&l2_cache_size, hipDevAttrL2CacheSize, device);

  int l2_block_size = (int) sqrt(0.75 * ((double) l2_cache_size / 8));
  int min_block_size = LDA / max_submat_blocks;

  int block_size = std::max(l2_block_size, min_block_size);
  block_size = std::min(block_size, LDA);

  return block_size;
}

uint32_t get_device_sm_count(int32_t device) {
  int num_sm;
  hipDeviceGetAttribute(&num_sm, hipDevAttrMultiProcessorCount, device);

  return num_sm;
}

}
}
