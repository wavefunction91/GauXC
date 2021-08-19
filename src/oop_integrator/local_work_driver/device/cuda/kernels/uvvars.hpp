#pragma once
#include "device/xc_device_task.hpp"

namespace GauXC {

void eval_uvvars_lda_cuda( size_t ntasks, int32_t nbe_max, int32_t npts_max,
  XCDeviceTask* device_tasks, cudaStream_t stream );

void eval_uvvars_gga_cuda( size_t ntasks, size_t npts_total, int32_t nbe_max, 
  int32_t npts_max, XCDeviceTask* device_tasks, const double* denx, 
  const double* deny, const double* denz, double* gamma, cudaStream_t stream );

}
