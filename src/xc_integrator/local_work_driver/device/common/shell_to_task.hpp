#pragma once
#include <gauxc/shell.hpp>

namespace GauXC {

struct ShellToTaskDevice {
  int32_t  ntask;
  int32_t*  task_idx_device;
  int32_t*  task_shell_offs_device;
  Shell<double>*   shell_device;
};

struct AngularMomentumShellToTaskBatch {
  size_t             ntask_average;
  size_t             npts_average;
  size_t             nshells_in_batch;
  ShellToTaskDevice* shell_to_task_device;
  uint32_t           pure;
};

}
