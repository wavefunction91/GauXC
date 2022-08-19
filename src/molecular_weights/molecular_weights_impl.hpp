#pragma once
#include <gauxc/molecular_weights.hpp>
#include <gauxc/xc_integrator/local_work_driver.hpp>

namespace GauXC::detail {

class MolecularWeightsImpl {

protected:

  std::unique_ptr<LocalWorkDriver> local_work_driver_;
  util::Timer timer_;

public:

  MolecularWeightsImpl() = delete;
  virtual ~MolecularWeightsImpl() noexcept = default;
  MolecularWeightsImpl( const MolecularWeightsImpl& ) = delete;
  MolecularWeightsImpl( MolecularWeightsImpl&& ) noexcept = default;

  inline MolecularWeightsImpl(std::unique_ptr<LocalWorkDriver>&& lwd) :
    local_work_driver_(std::move(lwd)) {}

  virtual void modify_weights(LoadBalancer&) const = 0;
  inline const util::Timer& get_timings() const {
    return timer_;
  };

  inline util::Timer& get_timer() {
    return timer_;
  };
};

}
