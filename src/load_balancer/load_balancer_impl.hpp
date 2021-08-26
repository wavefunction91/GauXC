#pragma once

#include <gauxc/load_balancer.hpp>

namespace GauXC  {
namespace detail {

class LoadBalancerImpl {

public:

  using basis_type = BasisSet<double>;

protected:

  GAUXC_MPI_CODE(MPI_Comm comm_;)

  std::shared_ptr<Molecule>   mol_;
  std::shared_ptr<MolGrid>    mg_;
  std::shared_ptr<basis_type> basis_;
  std::shared_ptr<MolMeta>    molmeta_;

  std::vector< XCTask >     local_tasks_;

  LoadBalancerState         state_;

  util::Timer               timer_;

  virtual std::vector< XCTask > create_local_tasks_() const = 0;

public:

  LoadBalancerImpl() = delete;

  LoadBalancerImpl( GAUXC_MPI_CODE(MPI_Comm,) const Molecule&, const MolGrid& mg,  
    const basis_type& );
  LoadBalancerImpl( GAUXC_MPI_CODE(MPI_Comm,) const Molecule&, const MolGrid& mg,  
    const basis_type&, const MolMeta& );
  LoadBalancerImpl( GAUXC_MPI_CODE(MPI_Comm,) const Molecule&, const MolGrid& mg,  
    const basis_type&, std::shared_ptr<MolMeta> );

  LoadBalancerImpl( const LoadBalancerImpl& );
  LoadBalancerImpl( LoadBalancerImpl&& ) noexcept;

  virtual ~LoadBalancerImpl() noexcept;

  const std::vector< XCTask >& get_tasks() const;
        std::vector< XCTask >& get_tasks()      ;

  const util::Timer& get_timings() const;

  size_t max_npts()       const;
  size_t max_nbe()        const;
  size_t max_npts_x_nbe() const;

  const Molecule& molecule() const;
  const MolMeta&  molmeta()  const;
  const basis_type& basis()  const;

  LoadBalancerState& state();

#ifdef GAUXC_ENABLE_MPI
  MPI_Comm comm() const;
#endif

  virtual std::unique_ptr<LoadBalancerImpl> clone() const = 0;

};

}
}
