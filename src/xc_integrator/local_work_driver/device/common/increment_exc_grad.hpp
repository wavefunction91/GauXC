/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include <gauxc/shell.hpp>
#include "device/xc_device_task.hpp"
#include "device/device_queue.hpp"
#include "shell_to_task.hpp"

namespace GauXC {

void increment_exc_grad_lda( size_t nshell, ShellToTaskDevice* shell_to_task,
  XCDeviceTask* device_tasks, double* EXC_GRAD, device_queue );
void increment_exc_grad_gga( size_t nshell, ShellToTaskDevice* shell_to_task,
  XCDeviceTask* device_tasks, double* EXC_GRAD, device_queue );

}

