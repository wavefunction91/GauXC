#pragma once
#include <gauxc/xc_task.hpp>

namespace GauXC      {
namespace integrator {
namespace hip       {

using namespace GauXC::hip;

template <typename T>
void zmat_lda_hip( size_t           ntasks,
                    int32_t          max_nbf,
                    int32_t          max_npts,
                    XCTaskDevice<T>* tasks_device,
                    hipStream_t     stream );

template <typename T>
void zmat_gga_hip( size_t           ntasks,
                    int32_t          max_nbf,
                    int32_t          max_npts,
                    XCTaskDevice<T>* tasks_device,
                    hipStream_t     stream );
              
}
}
}
