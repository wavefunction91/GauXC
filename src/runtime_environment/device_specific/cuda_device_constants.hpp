#pragma once
#include <cstdint>

namespace GauXC {
namespace cuda  {

static constexpr uint32_t warp_size = 32;
static constexpr uint32_t max_threads_per_thread_block = 1024;
static constexpr uint32_t max_warps_per_thread_block = 
  max_threads_per_thread_block / warp_size;

namespace obara_saika {

static constexpr int points_per_subtask = 256;
static constexpr int max_primpairs = 32;

static constexpr bool k0_use_shared = true;
static constexpr bool k1_use_shared = true;
static constexpr bool k2_use_shared = true;

static constexpr bool k00_use_shared = true;
static constexpr bool k10_use_shared = true;
static constexpr bool k11_use_shared = true;
static constexpr bool k20_use_shared = true;
static constexpr bool k21_use_shared = true;
static constexpr bool k22_use_shared = true;

}

}

}
