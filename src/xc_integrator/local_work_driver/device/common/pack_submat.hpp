#include "device/xc_device_task.hpp"
#include "device/device_queue.hpp"

namespace GauXC {

void sym_pack_submat( size_t ntasks, XCDeviceTask* device_tasks, const double* A,
  int32_t LDA, int32_t submat_block_size, device_queue queue );

void asym_pack_submat( size_t ntasks, XCDeviceTask* device_tasks, const double* A,
  int32_t LDA, int32_t submat_block_size, device_queue queue );
}
