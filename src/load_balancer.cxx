#include "load_balancer_impl.hpp"
#include "load_balancer_defaults.hpp"

namespace GauXC {

LoadBalancer::LoadBalancer( std::unique_ptr<pimpl_type>&& pimpl ): 
  pimpl_( std::move(pimpl) ) { }

LoadBalancer::LoadBalancer() : LoadBalancer( nullptr ) { }

LoadBalancer::LoadBalancer( MPI_Comm comm, const Molecule& mol, const MolGrid& mg, 
  const basis_type& basis, std::shared_ptr<MolMeta> meta ) : 
  LoadBalancer( detail::make_default_load_balancer( comm, mol, mg, basis, meta ) ) 
  { }

LoadBalancer::LoadBalancer( MPI_Comm comm, const Molecule& mol, const MolGrid& mg, 
  const basis_type& basis, const MolMeta& meta ) : 
  LoadBalancer( detail::make_default_load_balancer( comm, mol, mg, basis, meta ) ) 
  { }

LoadBalancer::LoadBalancer( MPI_Comm comm, const Molecule& mol, const MolGrid& mg, 
  const basis_type& basis ) : 
  LoadBalancer( detail::make_default_load_balancer( comm, mol, mg, basis ) ) 
  { }


LoadBalancer::LoadBalancer( const LoadBalancer& other ) :
  LoadBalancer(other.pimpl_->clone()){ };

LoadBalancer::LoadBalancer( LoadBalancer&& ) noexcept = default;
              
LoadBalancer::~LoadBalancer() noexcept = default;;


const std::vector<XCTask>& LoadBalancer::get_tasks() const {
  if( not pimpl_ ) throw std::runtime_error("Not Initialized");
  return pimpl_->get_tasks();
}
std::vector<XCTask>& LoadBalancer::get_tasks() {
  if( not pimpl_ ) throw std::runtime_error("Not Initialized");
  return pimpl_->get_tasks();
}

}
