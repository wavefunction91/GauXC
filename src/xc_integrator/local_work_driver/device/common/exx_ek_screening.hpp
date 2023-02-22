/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for retails
 */
#include "device/xc_device_task.hpp"
#include "device/device_queue.hpp"
#include <gauxc/shell.hpp>

namespace GauXC {

void exx_ek_screening_bfn_stats( size_t        ntasks,
                                 XCDeviceTask* tasks_device,
                                 double      * bfn_max_device,
                                 size_t        LDBFM,
                                 device_queue queue );

void exx_ek_collapse_fmax_to_shells(
  int                  ntask,
  int                  nshells,
  const Shell<double>* shells_device,
  const int32_t*       shell_to_bf,
  const double*        fmax_bfn_device,
  size_t               LDF_bfn,
  double*              fmax_shell_device,
  size_t               LDF_shell,
  device_queue         queue
);

}
