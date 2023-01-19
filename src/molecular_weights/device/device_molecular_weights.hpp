#pragma once
#include "../molecular_weights_impl.hpp"
namespace GauXC::detail {

class DeviceMolecularWeights : public MolecularWeightsImpl {

public:

  DeviceMolecularWeights() = delete;
  virtual ~DeviceMolecularWeights() noexcept = default;
  DeviceMolecularWeights( const DeviceMolecularWeights& ) = delete;
  DeviceMolecularWeights( DeviceMolecularWeights&& ) noexcept = default;

  inline DeviceMolecularWeights(std::unique_ptr<LocalWorkDriver>&& lwd) :
    MolecularWeightsImpl(std::move(lwd)) {}

  void modify_weights(LoadBalancer&) const final;

};

template <typename... Args>
std::unique_ptr<MolecularWeightsImpl> 
  make_device_mol_weights_impl(Args&&... args) {
  return std::make_unique<DeviceMolecularWeights>(
    std::forward<Args>(args)...);
}

}
