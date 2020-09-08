#pragma once

#include <SYCL/sycl.hpp>
#include <gauxc/xc_task.hpp>

namespace GauXC      {
namespace integrator {
namespace sycl       {

    using namespace GauXC::sycl;

    template <typename T>
    void task_pack_density_matrix(size_t ntasks, XCTaskDevice<T> *device_tasks,
                                  T *P_device, size_t LDP, cl::sycl::queue *stream);
}
}
}
