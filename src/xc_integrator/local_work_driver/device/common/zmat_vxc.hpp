/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for retails
 */
#include "device/xc_device_task.hpp"
#include "device/device_queue.hpp"

namespace GauXC {

void zmat_lda_vxc( size_t        ntasks,
                   int32_t       max_nbf,
                   int32_t       max_npts,
                   XCDeviceTask* tasks_device,
                   device_queue queue );

void zmat_gga_vxc( size_t        ntasks,
                   int32_t       max_nbf,
                   int32_t       max_npts,
                   XCDeviceTask* tasks_device,
                   device_queue queue );

}
