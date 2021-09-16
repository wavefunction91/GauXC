#include "device/xc_device_task.hpp"
#include "device/type_erased_queue.hpp"

namespace GauXC {

void pack_submat( size_t ntasks, XCDeviceTask* device_tasks, const double* A,
  int32_t LDA, int32_t submat_block_size, type_erased_queue queue );

}
