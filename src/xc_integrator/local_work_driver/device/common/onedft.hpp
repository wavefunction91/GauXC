/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include "device/xc_device_task.hpp"
#include "device/xc_device_data.hpp"
#include "device/device_queue.hpp"

namespace GauXC {

void sz_to_ab(  size_t sz,
                const void* src_a,
                const void* src_b,
                void* dest_a,
                void* dest_b,
                device_queue queue );

void zmat_onedft_vxc( size_t          ntasks,
    int32_t           max_nbf,
    int32_t           max_npts,
    XCDeviceTask*     tasks_device,
    integrator_xc_approx scheme,
    density_id sel,
    device_queue queue );
}
