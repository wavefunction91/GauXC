#pragma once
#include <cstdint>

namespace GauXC {
namespace cuda  {

static constexpr uint32_t warp_size = 32;
static constexpr uint32_t max_threads_per_thread_block = 1024;
static constexpr uint32_t max_warps_per_thread_block = 
  max_threads_per_thread_block / warp_size;

}
}
