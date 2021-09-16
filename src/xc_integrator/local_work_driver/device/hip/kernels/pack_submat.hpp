#pragma once
#include <hip/hip_runtime.h>
#include "device/xc_device_task.hpp"

namespace GauXC {

void pack_submat( size_t ntasks, XCDeviceTask* device_tasks, const double* A,
  int32_t LDA, int32_t submat_block_size, hipStream_t stream );

}
