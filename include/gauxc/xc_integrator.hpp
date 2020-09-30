#pragma once

#include <memory>

#include <gauxc/types.hpp>
#include <gauxc/load_balancer.hpp>

namespace GauXC {

namespace detail {
  template <typename MatrixType>
  class XCIntegratorImpl;
}


template <typename MatrixType>
class XCIntegrator {

public:

  using matrix_type   = MatrixType;
  using value_type    = typename matrix_type::value_type;  
  using basisset_type = BasisSet< value_type >;

  using exc_vxc_type = std::tuple< value_type, matrix_type >;

private:

  using pimpl_type    = detail::XCIntegratorImpl<MatrixType>;

  std::unique_ptr<pimpl_type> pimpl_;

public:

  XCIntegrator();
  XCIntegrator( std::unique_ptr<pimpl_type>&& pimpl );

#ifdef GAUXC_ENABLE_MPI

  XCIntegrator( ExecutionSpace, MPI_Comm, const functional_type&, 
                const basisset_type&, std::shared_ptr<LoadBalancer> );
  XCIntegrator( MPI_Comm, const functional_type&, const basisset_type&, 
                std::shared_ptr<LoadBalancer> );
  //XCIntegrator( MPI_Comm, const functional_type&, const basisset_type&,
  //  LoadBalancer&& );

  //XCIntegrator( MPI_Comm, const functional_type&, const basisset_type&, 
  //  const Molecule&, const MolGrid& );
  //XCIntegrator( MPI_Comm, const functional_type&, const basisset_type&, 
  //  const Molecule&, const MolGrid&, const MolMeta );
  //XCIntegrator( MPI_Comm, const functional_type&, const basisset_type&, 
  //  const Molecule&, const MolGrid&, const std::shared_ptr<MolMeta> );

#else

  XCIntegrator( ExecutionSpace, const functional_type&, 
                const basisset_type&, std::shared_ptr<LoadBalancer> );
  XCIntegrator( const functional_type&, const basisset_type&, 
                std::shared_ptr<LoadBalancer> );
  //XCIntegrator( const functional_type&, const basisset_type&,
  //  LoadBalancer&& );

  //XCIntegrator( const functional_type&, const basisset_type&, 
  //  const Molecule&, const MolGrid& );
  //XCIntegrator( const functional_type&, const basisset_type&, 
  //  const Molecule&, const MolGrid&, const MolMeta );
  //XCIntegrator( const functional_type&, const basisset_type&, 
  //  const Molecule&, const MolGrid&, const std::shared_ptr<MolMeta> );

#endif
  ~XCIntegrator() noexcept;


  exc_vxc_type eval_exc_vxc( const MatrixType& );
};


}
