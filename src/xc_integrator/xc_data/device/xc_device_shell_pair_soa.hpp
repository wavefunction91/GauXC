#pragma once
#include <gauxc/shell_pair.hpp>

namespace GauXC {
struct XCDeviceShellPairSoA {
  using shell_pair = ShellPair<double>;
  using point      = detail::cartesian_point;
  std::vector<shell_pair*>            shell_pair_dev_ptr;
  std::vector<int32_t>                shell_pair_nprim_pairs;
  std::vector<std::pair<int,int>>     shell_pair_shidx;
  std::vector<std::pair<int,int>>     shell_pair_ls;
  std::vector<std::pair<point,point>> shell_pair_centers;

  inline void reset() {
    shell_pair_nprim_pairs.clear();
    shell_pair_dev_ptr.clear();
    shell_pair_ls.clear();
    shell_pair_centers.clear();
  }
};
}
