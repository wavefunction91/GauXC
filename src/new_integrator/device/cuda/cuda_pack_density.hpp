#pragma once
#include <gauxc/xc_task.hpp>

namespace GauXC      {
namespace integrator {
namespace cuda       {

using namespace GauXC::cuda;

template <typename T>
void task_pack_density_matrix( size_t           ntasks,
                               XCTaskDevice<T>* device_tasks,
                               T*               P_device,
                               size_t           LDP,
                               cudaStream_t     stream );
                               
}
}
}
