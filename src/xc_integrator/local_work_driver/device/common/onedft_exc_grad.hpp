/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy).
 *
 * (c) 2024-2025, Microsoft Corporation
 *
 * All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include "device/xc_device_task.hpp"
#include "device/device_queue.hpp"

namespace GauXC {

// Transform OneDFT per-direction Vxc (stored in gamma_pp/vgamma_pp etc.)
// into the standard vrho/vgamma/dden format expected by inc_exc_grad kernels.
// After this call, dden_sx/sy/sz contain vds_x/y/z and vgamma_pp/pm/mm = 1/0/1.
void transform_onedft_vxc_for_grad(
    size_t ntasks,
    int32_t max_npts,
    XCDeviceTask* tasks_device,
    device_queue queue );

}
