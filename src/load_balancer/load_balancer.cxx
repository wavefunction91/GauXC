#include "load_balancer_impl.hpp"

namespace GauXC {

LoadBalancer::LoadBalancer( std::unique_ptr<pimpl_type>&& pimpl ): 
  pimpl_( std::move(pimpl) ) { }

LoadBalancer::LoadBalancer() : LoadBalancer( nullptr ) { }

LoadBalancer::LoadBalancer( const LoadBalancer& other ) :
  LoadBalancer(other.pimpl_->clone()){ }

LoadBalancer::LoadBalancer( LoadBalancer&& ) noexcept = default;
              
LoadBalancer::~LoadBalancer() noexcept = default;


const std::vector<XCTask>& LoadBalancer::get_tasks() const {
  if( not pimpl_ ) throw std::runtime_error("Not Initialized");
  return pimpl_->get_tasks();
}
std::vector<XCTask>& LoadBalancer::get_tasks() {
  if( not pimpl_ ) throw std::runtime_error("Not Initialized");
  return pimpl_->get_tasks();
}

const util::Timer& LoadBalancer::get_timings() const {
  if( not pimpl_ ) throw std::runtime_error("Not Initialized");
  return pimpl_->get_timings();
}

size_t LoadBalancer::max_npts() const {
  if( not pimpl_ ) throw std::runtime_error("Not Initialized");
  return pimpl_->max_npts();
}
size_t LoadBalancer::max_nbe() const {
  if( not pimpl_ ) throw std::runtime_error("Not Initialized");
  return pimpl_->max_nbe();
}
size_t LoadBalancer::max_npts_x_nbe() const {
  if( not pimpl_ ) throw std::runtime_error("Not Initialized");
  return pimpl_->max_npts_x_nbe();
}



const Molecule& LoadBalancer::molecule() const {
  if( not pimpl_ ) throw std::runtime_error("Not Initialized");
  return pimpl_->molecule();
}
const MolMeta& LoadBalancer::molmeta() const {
  if( not pimpl_ ) throw std::runtime_error("Not Initialized");
  return pimpl_->molmeta();
}

const LoadBalancer::basis_type& LoadBalancer::basis() const {
  if( not pimpl_ ) throw std::runtime_error("Not Initialized");
  return pimpl_->basis();
}

LoadBalancerState& LoadBalancer::state() {
  if( not pimpl_ ) throw std::runtime_error("Not Initialized");
  return pimpl_->state();
}

#ifdef GAUXC_ENABLE_MPI
MPI_Comm LoadBalancer::comm() const {
  if( not pimpl_ ) throw std::runtime_error("Not Initialized");
  return pimpl_->comm();
}
#endif

}
