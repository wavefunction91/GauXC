#pragma once
#include <gauxc/xc_task.hpp>
#include <host/local_host_work_driver.hpp>

namespace GauXC {
namespace exx_detail {
  using host_task_container = std::vector<XCTask>;
  using host_task_iterator  = typename host_task_container::iterator;
}

std::vector<std::vector<int32_t>> exx_ek_screening( 
  const BasisSet<double>& basis, const BasisSetMap& basis_map,
  const double* P_abs, size_t ldp, const double* V_shell_max, size_t ldv,
  double eps_E, double eps_K, LocalHostWorkDriver* lwd, 
  exx_detail::host_task_iterator task_begin,
  exx_detail::host_task_iterator task_end );

}
