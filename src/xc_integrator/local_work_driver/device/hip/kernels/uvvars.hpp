#pragma once
#include "hip/hip_runtime.h"
#include "device/xc_device_task.hpp"

namespace GauXC {

void eval_uvvars_lda_hip( size_t ntasks, int32_t nbe_max, int32_t npts_max,
  XCDeviceTask* device_tasks, hipStream_t stream );

void eval_uvvars_gga_hip( size_t ntasks, size_t npts_total, int32_t nbe_max, 
  int32_t npts_max, XCDeviceTask* device_tasks, const double* denx, 
  const double* deny, const double* denz, double* gamma, hipStream_t stream );

}
