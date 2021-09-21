#pragma once
#include <gauxc/shell.hpp>

namespace GauXC {

struct ShellToTaskDevice {
  int32_t  ntask;
  int32_t*  task_idx_device;
  int32_t*  task_shell_offs_device;
  Shell<double>*   shell_device;
};

}
