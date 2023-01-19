#pragma once
#include "device/xc_device_task.hpp"
#include "device/device_queue.hpp"

namespace GauXC {

void sym_task_inc_potential( size_t        ntasks,
                         XCDeviceTask* device_tasks,
                         double*       V_device,
                         size_t        LDV,
                         size_t        submat_block,
                         device_queue  queue );
                               
void asym_task_inc_potential( size_t        ntasks,
                         XCDeviceTask* device_tasks,
                         double*       V_device,
                         size_t        LDV,
                         size_t        submat_block,
                         device_queue  queue );
}

