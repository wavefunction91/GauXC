/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include "device/xc_device_task.hpp"
#include "device/xc_device_data.hpp"
#include "device/device_queue.hpp"

namespace GauXC {


void eval_uvars_lda( size_t ntasks, int32_t npts_max, integrator_ks_scheme ks_scheme,
  XCDeviceTask* device_tasks, device_queue queue );

void eval_uvars_gga( size_t ntasks, int32_t npts_max, integrator_ks_scheme ks_scheme,
  XCDeviceTask* device_tasks, device_queue queue );

void eval_vvar( size_t ntasks, int32_t nbf_max, int32_t npts_max, bool do_grad, density_id den_select,
  XCDeviceTask* device_tasks, device_queue queue );
}
