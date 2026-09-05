#pragma once
#include <gauxc/xc_task.hpp>

namespace GauXC      {
namespace integrator {
namespace hip       {

using namespace GauXC::hip;

template <typename T>
void task_pack_density_matrix( size_t           ntasks,
                               XCTaskDevice<T>* device_tasks,
                               T*               P_device,
                               size_t           LDP,
                               hipStream_t     stream );
                               
}
}
}
