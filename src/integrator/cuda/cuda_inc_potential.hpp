#pragma once
#include <gauxc/xc_task.hpp>

namespace GauXC      {
namespace integrator {
namespace cuda       {

using namespace GauXC::cuda;

template <typename T>
void task_inc_potential( size_t           ntasks,
                         XCTaskDevice<T>* device_tasks,
                         T*               V_device,
                         size_t           LDV,
                         cudaStream_t     stream );

template <typename T>
void symmetrize_matrix( size_t       nbf,
                        size_t       LDV,
                        T*           V_device,
                        cudaStream_t stream);

}
}
}

