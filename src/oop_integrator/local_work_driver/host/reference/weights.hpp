#include "host/local_host_work_driver_pimpl.hpp"

namespace GauXC {

using task_iterator = detail::LocalHostWorkDriverPIMPL::task_iterator;

void reference_ssf_weights_host(
  const Molecule&        mol,
  const MolMeta&         meta,
  task_iterator          task_begin,
  task_iterator          task_end
);

void reference_becke_weights_host(
  const Molecule&        mol,
  const MolMeta&         meta,
  task_iterator          task_begin,
  task_iterator          task_end
);

}
