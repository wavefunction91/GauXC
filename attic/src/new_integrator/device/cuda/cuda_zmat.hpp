#pragma once
#include <gauxc/xc_task.hpp>

namespace GauXC      {
namespace integrator {
namespace cuda       {

using namespace GauXC::cuda;

template <typename T>
void zmat_lda_cuda( size_t           ntasks,
                    int32_t          max_nbf,
                    int32_t          max_npts,
                    XCTaskDevice<T>* tasks_device,
                    cudaStream_t     stream );

template <typename T>
void zmat_gga_cuda( size_t           ntasks,
                    int32_t          max_nbf,
                    int32_t          max_npts,
                    XCTaskDevice<T>* tasks_device,
                    cudaStream_t     stream );
              
}
}
}
