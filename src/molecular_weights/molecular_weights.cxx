#include <gauxc/exceptions.hpp>
#include "molecular_weights_impl.hpp"
#include "host/host_molecular_weights.hpp"
#ifdef GAUXC_ENABLE_DEVICE
#include "device/device_molecular_weights.hpp"
#endif

namespace GauXC {

MolecularWeights::~MolecularWeights() noexcept = default;
MolecularWeights::MolecularWeights(MolecularWeights&&) noexcept = default;

MolecularWeights::MolecularWeights(pimpl_ptr_type&& ptr) :
  pimpl_(std::move(ptr)) {}

void MolecularWeights::modify_weights(load_balancer_reference lb) const {
  if(not pimpl_) GAUXC_PIMPL_NOT_INITIALIZED();
  pimpl_->modify_weights(lb);
}

const util::Timer& MolecularWeights::get_timings() const {
  if(not pimpl_) GAUXC_PIMPL_NOT_INITIALIZED();
  return pimpl_->get_timings();
}







MolecularWeightsFactory::MolecularWeightsFactory( ExecutionSpace ex,
  std::string lwd_kernel, MolecularWeightsSettings settings ) :
  ex_(ex), lwd_kernel_(lwd_kernel), settings_(settings) {}

std::shared_ptr<MolecularWeights> 
  MolecularWeightsFactory::get_shared_instance() {

  // Create Local Work Driver
  LocalWorkSettings lwd_settings;
  auto lwd = LocalWorkDriverFactory::make_local_work_driver(ex_,
    lwd_kernel_, lwd_settings );

  // Create MolecularWeights instance
  if( ex_ == ExecutionSpace::Host ) {
    return std::make_shared<MolecularWeights>(
      detail::make_host_mol_weights_impl(std::move(lwd))
    );
  } else {
  #ifdef GAUXC_ENABLE_DEVICE
    return std::make_shared<MolecularWeights>(
      detail::make_device_mol_weights_impl(std::move(lwd))
    );
  #else
    GAUXC_GENERIC_EXCEPTION("Device API Not Enabled");
  #endif
  }

}

}
