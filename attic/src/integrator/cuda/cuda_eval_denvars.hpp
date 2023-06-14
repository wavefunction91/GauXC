#pragma once
#include <gauxc/xc_task.hpp>

namespace GauXC      {
namespace integrator {
namespace cuda       {

using namespace GauXC::cuda;

template <typename T>
void eval_uvars_lda_device( size_t           ntasks,
                            size_t           max_nbf,
                            size_t           max_npts,
                            XCTaskDevice<T>* tasks_device,
                            cudaStream_t     stream );

template <typename T>
void eval_uvars_gga_device( size_t           ntasks,
                            size_t           max_nbf,
                            size_t           max_npts,
                            XCTaskDevice<T>* tasks_device,
                            cudaStream_t     stream );
 

template <typename T>
void eval_vvars_gga_device( size_t       npts,
                            const T*     den_x_device,
                            const T*     den_y_device,
                            const T*     den_z_device,
                                  T*     gamma_device,
                            cudaStream_t stream );
                          

}
}
}
