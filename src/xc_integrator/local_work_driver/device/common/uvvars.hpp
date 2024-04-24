/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include "device/xc_device_task.hpp"
#include "device/device_queue.hpp"

namespace GauXC {

void eval_uvvars_lda( size_t ntasks, int32_t nbe_max, int32_t npts_max,
  XCDeviceTask* device_tasks, device_queue queue );

void eval_uvvars_gga( size_t ntasks, size_t npts_total, int32_t nbe_max, 
  int32_t npts_max, XCDeviceTask* device_tasks, const double* denx, 
  const double* deny, const double* denz, double* gamma, device_queue queue );

void eval_uvvars_mgga( size_t ntasks, size_t npts_total, int32_t nbe_max, 
  int32_t npts_max, XCDeviceTask* device_tasks, const double* denx, 
  const double* deny, const double* denz, double* gamma, bool do_lapl,
  device_queue queue );

}
