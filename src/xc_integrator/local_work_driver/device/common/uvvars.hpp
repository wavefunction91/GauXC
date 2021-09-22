#pragma once
#include "device/xc_device_task.hpp"
#include "device/device_queue.hpp"

namespace GauXC {

void eval_uvvars_lda( size_t ntasks, int32_t nbe_max, int32_t npts_max,
  XCDeviceTask* device_tasks, device_queue queue );

void eval_uvvars_gga( size_t ntasks, size_t npts_total, int32_t nbe_max, 
  int32_t npts_max, XCDeviceTask* device_tasks, const double* denx, 
  const double* deny, const double* denz, double* gamma, device_queue queue );

}
