#pragma once
#include <cstdint>

namespace GauXC {
namespace cuda  {

static constexpr uint32_t warp_size = 32;
static constexpr uint32_t max_threads_per_thread_block = 1024;
static constexpr uint32_t max_warps_per_thread_block = 
  max_threads_per_thread_block / warp_size;

static constexpr uint32_t max_submat_blocks = 10;

// Properties for weight algorithm
static constexpr uint32_t weight_unroll = 2;
static_assert(weight_unroll == 2, "Weight unroll is only tested for value of 2");
static constexpr uint32_t weight_thread_block = 640;
static constexpr uint32_t weight_thread_block_per_sm = 2;

uint32_t get_submat_cut_block(int32_t LDA);
uint32_t get_device_sm_count();
}
}
