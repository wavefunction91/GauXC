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
                                 double      * max_bfn_sum_device,
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

void exx_ek_shellpair_collision(
  int32_t       ntasks,
  int32_t       nshells,
  const double* V_max_device,
  size_t        LDV,
  const double* F_max_shl_device,
  size_t        LDF,
  const double* max_bf_sum_device,
  double        eps_E,
  double        eps_K,
  uint32_t*     collisions,
  int           LD_coll,
  uint32_t*     counts,
  void*         dyn_stack,
  size_t        dyn_size,
  device_queue  queue
);

}
