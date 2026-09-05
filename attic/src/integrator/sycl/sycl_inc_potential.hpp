#pragma once

#include <CL/sycl.hpp>
#include <gauxc/xc_task.hpp>

namespace GauXC      {
namespace integrator {
namespace sycl       {

using namespace GauXC::sycl;

template <typename T>
void task_inc_potential(size_t ntasks, XCTaskDevice<T> *device_tasks,
                        T *V_device, size_t LDV, cl::sycl::queue *queue);
}
}
}

