#pragma once 
#include <gauxc/shell_pair.hpp>

namespace GauXC {

struct ShellPairToTaskHost {
  using shell_pair = ShellPair<double>;
  using point      = detail::cartesian_point;
  std::vector<int32_t> task_idx;
  std::vector<int32_t> task_shell_off_row;
  std::vector<int32_t> task_shell_off_col;
  shell_pair* shell_pair_device;

  int32_t idx_bra, idx_ket;
  int32_t lA, lB;
  point rA, rB;

  void clear() {
    task_idx.clear();
    task_shell_off_row.clear();
    task_shell_off_col.clear();
    shell_pair_device = nullptr;
  }
};

struct ShellPairToTaskDevice {
  using shell_pair = ShellPair<double>;
  using point      = detail::cartesian_point;
  int32_t* task_idx_device;
  int32_t* task_shell_off_row_device;
  int32_t* task_shell_off_col_device;
  shell_pair* shell_pair_device;
  int32_t ntask;

  double X_AB, Y_AB, Z_AB;
};

struct AngularMomentumShellPairToTaskBatch {
  size_t ntask_average;
  size_t npts_average;
  size_t nshells_in_batch;
  ShellPairToTaskDevice* shell_pair_to_task_device;

  int lA, lB;
};

}
