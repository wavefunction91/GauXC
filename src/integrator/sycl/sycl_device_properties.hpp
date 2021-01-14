#pragma once
#include <CL/sycl.hpp>

namespace GauXC {
namespace sycl  {

    //static cl::sycl::device dev(cl::sycl::gpu_selector{});
    // cl::sycl::queue q(dev);
    // cl::sycl::buffer<int, 1> buf(cl::sycl::range<1>(1));
    // cl::sycl::program prg(q.get_context());
    // prg.build_with_kernel_type<class SingleTask>();
    // assert(prg.has_kernel<class SingleTask>());
    // cl::sycl::kernel krn = prg.get_kernel<class SingleTask>();
    //
    // q.submit([&](cl::sycl::handler &cgh) {
    //     auto acc = buf.get_access<cl::sycl::access::mode::read_write>(cgh);
    //     cgh.single_task<class SingleTask>(krn, [=]() { acc[0] = acc[0] + 1; }); });

    static size_t warp_size = 32; //krn.get_info<cl::sycl::info::kernel_work_group::preferred_work_group_size_multiple>(dev);
    static size_t max_threads_per_thread_block = 256;
    static size_t max_warps_per_thread_block = max_threads_per_thread_block / warp_size;

}
}
