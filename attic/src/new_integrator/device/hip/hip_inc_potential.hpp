#pragma once
#include <gauxc/xc_task.hpp>

namespace GauXC      {
namespace integrator {
namespace hip       {

using namespace GauXC::hip;

template <typename T>
void task_inc_potential( size_t           ntasks,
                         XCTaskDevice<T>* device_tasks,
                         T*               V_device,
                         size_t           LDV,
                         hipStream_t     stream );
                               
}
}
}

