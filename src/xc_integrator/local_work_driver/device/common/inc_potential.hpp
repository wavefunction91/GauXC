#pragma once
#include "device/xc_device_task.hpp"
#include "device/type_erased_queue.hpp"

namespace GauXC {

void task_inc_potential( size_t        ntasks,
                         XCDeviceTask* device_tasks,
                         double*       V_device,
                         size_t        LDV,
                         size_t        submat_block,
                         type_erased_queue  queue );
                               
}

