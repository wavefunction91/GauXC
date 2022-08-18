#include "host_molecular_weights.hpp"
#include "host/local_host_work_driver.hpp"

namespace GauXC::detail {

void HostMolecularWeights::modify_weights( LoadBalancer& lb ) const {

  if(lb.state().modified_weights_are_stored)
    GAUXC_GENERIC_EXCEPTION("Attempting to Overwrite Modified Weights");

  // Cast LWD to LocalHostWorkDriver
  auto* lwd = dynamic_cast<LocalHostWorkDriver*>(this->local_work_driver_.get());

  // (Possibly) Generate tasks
  auto& tasks = lb.get_tasks();

  // Sort tasks on size (XXX: maybe doesnt matter?)
  auto task_comparator = []( const XCTask& a, const XCTask& b ) {
    return (a.points.size() * a.nbe) > (b.points.size() * b.nbe);
  };
  std::sort( tasks.begin(), tasks.end(), task_comparator );

  // Modify the weights
  const auto& mol  = lb.molecule();
  const auto& meta = lb.molmeta();
  lwd->partition_weights( XCWeightAlg::SSF, mol, meta, 
    tasks.begin(), tasks.end() );

  lb.state().modified_weights_are_stored = true;
}

}
