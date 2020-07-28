#pragma once

#include <mpi.h>
#include <gauxc/molgrid.hpp>
#include <gauxc/molmeta.hpp>
#include <gauxc/basisset.hpp>
#include <gauxc/xc_task.hpp>

namespace GauXC {

namespace detail {
  class LoadBalancerImpl;
}




class LoadBalancer {

  using basis_type = BasisSet<double>;
  using pimpl_type = detail::LoadBalancerImpl;
  std::unique_ptr<pimpl_type> pimpl_;

public:

  LoadBalancer();
  LoadBalancer( MPI_Comm, const Molecule&, const MolGrid&, const basis_type& );
  LoadBalancer( MPI_Comm, const Molecule&, const MolGrid&, const basis_type&,
    const MolMeta& );
  LoadBalancer( MPI_Comm, const Molecule&, const MolGrid&, const basis_type&,
    std::shared_ptr<MolMeta> );


  LoadBalancer( std::unique_ptr<pimpl_type>&& pimpl );

  LoadBalancer( const LoadBalancer& );
  LoadBalancer( LoadBalancer&& ) noexcept;

  ~LoadBalancer() noexcept;

  const std::vector<XCTask>& get_tasks() const;
        std::vector<XCTask>& get_tasks()      ;
  

  size_t max_npts()       const;
  size_t max_nbe()        const;
  size_t max_npts_x_nbe() const;

  const Molecule& molecule() const;
  const MolMeta&  molmeta()  const;

};





// Factories

namespace factory {

std::shared_ptr<LoadBalancer> make_default_load_balancer(
  MPI_Comm, const Molecule&, const MolGrid&, const BasisSet<double>&
);
std::shared_ptr<LoadBalancer> make_default_load_balancer(
  MPI_Comm, const Molecule&, const MolGrid&, const BasisSet<double>&, 
  const MolMeta&
);
std::shared_ptr<LoadBalancer> make_default_load_balancer(
  MPI_Comm, const Molecule&, const MolGrid&, const BasisSet<double>&, 
  std::shared_ptr<MolMeta>
);

}

}
