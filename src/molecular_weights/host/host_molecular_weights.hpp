#pragma once
#include "../molecular_weights_impl.hpp"
namespace GauXC::detail {

class HostMolecularWeights : public MolecularWeightsImpl {

public:

  HostMolecularWeights() = delete;
  virtual ~HostMolecularWeights() noexcept = default;
  HostMolecularWeights( const HostMolecularWeights& ) = delete;
  HostMolecularWeights( HostMolecularWeights&& ) noexcept = default;

  inline HostMolecularWeights(std::unique_ptr<LocalWorkDriver>&& lwd) :
    MolecularWeightsImpl(std::move(lwd)) {}

  void modify_weights(LoadBalancer&) const final;

};

template <typename... Args>
std::unique_ptr<MolecularWeightsImpl> 
  make_host_mol_weights_impl(Args&&... args) {
  return std::make_unique<HostMolecularWeights>(
    std::forward<Args>(args)...);
}

}
