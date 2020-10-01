#include "load_balancer_impl.hpp"
#include "load_balancer_defaults.hpp"

namespace GauXC {

LoadBalancer::LoadBalancer( std::unique_ptr<pimpl_type>&& pimpl ): 
  pimpl_( std::move(pimpl) ) { }

LoadBalancer::LoadBalancer() : LoadBalancer( nullptr ) { }

#ifdef GAUXC_ENABLE_MPI

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

#else

LoadBalancer::LoadBalancer( const Molecule& mol, const MolGrid& mg, 
  const basis_type& basis, std::shared_ptr<MolMeta> meta ) : 
  LoadBalancer( detail::make_default_load_balancer( mol, mg, basis, meta ) ) 
  { }

LoadBalancer::LoadBalancer( const Molecule& mol, const MolGrid& mg, 
  const basis_type& basis, const MolMeta& meta ) : 
  LoadBalancer( detail::make_default_load_balancer( mol, mg, basis, meta ) ) 
  { }

LoadBalancer::LoadBalancer( const Molecule& mol, const MolGrid& mg, 
  const basis_type& basis ) : 
  LoadBalancer( detail::make_default_load_balancer( mol, mg, basis ) ) 
  { }

#endif


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



namespace factory {

template <typename... Args>
std::shared_ptr<LoadBalancer> fwd_to_default( Args&&... args ) {
  return std::make_shared<LoadBalancer>( 
    detail::make_default_load_balancer( std::forward<Args>(args)... )
  );
}

#ifdef GAUXC_ENABLE_MPI

std::shared_ptr<LoadBalancer> make_default_load_balancer(
  MPI_Comm comm, const Molecule& mol, const MolGrid& mg, 
  const BasisSet<double> &basis
) {
  return fwd_to_default( comm, mol, mg, basis );
}


std::shared_ptr<LoadBalancer> make_default_load_balancer(
  MPI_Comm comm, const Molecule& mol, const MolGrid& mg, 
  const BasisSet<double> &basis, const MolMeta& meta
) {
  return fwd_to_default( comm, mol, mg, basis, meta );
}


std::shared_ptr<LoadBalancer> make_default_load_balancer(
  MPI_Comm comm, const Molecule& mol, const MolGrid& mg, 
  const BasisSet<double> &basis, std::shared_ptr<MolMeta> meta
) {
  return fwd_to_default( comm, mol, mg, basis, meta );
}

#else

std::shared_ptr<LoadBalancer> make_default_load_balancer(
  const Molecule& mol, const MolGrid& mg, 
  const BasisSet<double> &basis
) {
  return fwd_to_default( mol, mg, basis );
}


std::shared_ptr<LoadBalancer> make_default_load_balancer(
  const Molecule& mol, const MolGrid& mg, 
  const BasisSet<double> &basis, const MolMeta& meta
) {
  return fwd_to_default( mol, mg, basis, meta );
}


std::shared_ptr<LoadBalancer> make_default_load_balancer(
  const Molecule& mol, const MolGrid& mg, 
  const BasisSet<double> &basis, std::shared_ptr<MolMeta> meta
) {
  return fwd_to_default( mol, mg, basis, meta );
}

#endif


}

}
