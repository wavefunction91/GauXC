#include "load_balancer_impl.hpp"

namespace GauXC::detail {

#ifdef GAUXC_ENABLE_MPI

LoadBalancerImpl::LoadBalancerImpl( MPI_Comm comm, const Molecule& mol, 
  const MolGrid& mg, const basis_type& basis, std::shared_ptr<MolMeta> molmeta ) :
  comm_(comm), 
  mol_( std::make_shared<Molecule>(mol) ),
  mg_( std::make_shared<MolGrid>(mg)  ),
  basis_( std::make_shared<basis_type>(basis) ),
  molmeta_( molmeta ) { }

LoadBalancerImpl::LoadBalancerImpl( MPI_Comm comm, const Molecule& mol, 
  const MolGrid& mg, const basis_type& basis, const MolMeta& molmeta ) :
  LoadBalancerImpl( comm, mol, mg, basis, std::make_shared<MolMeta>(molmeta) ) { }

LoadBalancerImpl::LoadBalancerImpl( MPI_Comm comm, const Molecule& mol, 
  const MolGrid& mg, const basis_type& basis ) :
  LoadBalancerImpl( comm, mol, mg, basis, std::make_shared<MolMeta>(mol) ) { }

#else

LoadBalancerImpl::LoadBalancerImpl( const Molecule& mol, 
  const MolGrid& mg, const basis_type& basis, std::shared_ptr<MolMeta> molmeta ) :
  mol_( std::make_shared<Molecule>(mol) ),
  mg_( std::make_shared<MolGrid>(mg)  ),
  basis_( std::make_shared<basis_type>(basis) ),
  molmeta_( molmeta ) { }

LoadBalancerImpl::LoadBalancerImpl( const Molecule& mol, 
  const MolGrid& mg, const basis_type& basis, const MolMeta& molmeta ) :
  LoadBalancerImpl( mol, mg, basis, std::make_shared<MolMeta>(molmeta) ) { }

LoadBalancerImpl::LoadBalancerImpl( const Molecule& mol, 
  const MolGrid& mg, const basis_type& basis ) :
  LoadBalancerImpl( mol, mg, basis, std::make_shared<MolMeta>(mol) ) { }

#endif

LoadBalancerImpl::LoadBalancerImpl( const LoadBalancerImpl& ) = default;
LoadBalancerImpl::LoadBalancerImpl( LoadBalancerImpl&& ) noexcept = default;

LoadBalancerImpl::~LoadBalancerImpl() noexcept = default;

const std::vector<XCTask>& LoadBalancerImpl::get_tasks() const {
  if( not local_tasks_.size() ) throw std::runtime_error("No Tasks Created");
  return local_tasks_;
}

std::vector<XCTask>& LoadBalancerImpl::get_tasks() {
  auto create_tasks_st = std::chrono::high_resolution_clock::now();

  if( not local_tasks_.size() ) local_tasks_ = create_local_tasks_();

  auto create_tasks_en = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> create_tasks_dr = create_tasks_en - create_tasks_st; 
  timer_.add_timing("LoadBalancer.CreateTasks", create_tasks_dr);

  return local_tasks_;
}

const util::Timer& LoadBalancerImpl::get_timings() const {
  return timer_;
}


size_t LoadBalancerImpl::max_npts() const {

  if( not local_tasks_.size() ) return 0ul;

  return std::max_element( local_tasks_.cbegin(), local_tasks_.cend(),
    []( const auto& a, const auto& b ) {
      return a.points.size() < b.points.size();
    })->points.size();

}
size_t LoadBalancerImpl::max_nbe() const {

  if( not local_tasks_.size() ) return 0ul;

  return std::max_element( local_tasks_.cbegin(), local_tasks_.cend(),
    []( const auto& a, const auto& b ) {
      return a.nbe < b.nbe;
    })->nbe;

}
size_t LoadBalancerImpl::max_npts_x_nbe() const {

  if( not local_tasks_.size() ) return 0ul;

  auto it = std::max_element( local_tasks_.cbegin(), local_tasks_.cend(),
    []( const auto& a, const auto& b ) {
      return a.nbe * a.points.size() < b.nbe * b.points.size();
    });

  return it->nbe * it->points.size();

}



const Molecule& LoadBalancerImpl::molecule() const {
  return *mol_;
}

const MolMeta& LoadBalancerImpl::molmeta() const {
  return *molmeta_;
}

}
