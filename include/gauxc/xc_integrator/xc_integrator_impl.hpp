#pragma once

#include <gauxc/xc_integrator.hpp>
#include "xc_integrator_state.hpp"

namespace GauXC {
namespace detail {

template <typename MatrixType>
class XCIntegratorImpl {

public:

  using matrix_type   = MatrixType;
  using value_type    = typename matrix_type::value_type;
  using basisset_type = typename XCIntegrator<MatrixType>::basisset_type;
  using exc_vxc_type  = typename XCIntegrator<MatrixType>::exc_vxc_type;

protected:

  MPI_Comm comm_;
  std::shared_ptr<functional_type> func_;
  std::shared_ptr<basisset_type>   basis_;

  std::shared_ptr<LoadBalancer>    load_balancer_;
  XCIntegratorState                integrator_state_;

  virtual exc_vxc_type eval_exc_vxc_( const MatrixType& ) = 0;
  
public:

  XCIntegratorImpl( MPI_Comm                         comm, 
                    std::shared_ptr<functional_type> func, 
                    std::shared_ptr<basisset_type>   basis,
                    std::shared_ptr<LoadBalancer>    lb 
  ) : comm_(comm), func_(func), basis_(basis), load_balancer_(lb) { };

  template <typename... Args>
  XCIntegratorImpl( MPI_Comm                         comm, 
                    std::shared_ptr<functional_type> func, 
                    std::shared_ptr<basisset_type>   basis,
                    Args&&...                        args
  ) : comm_(comm), func_(func), basis_(basis), 
      load_balancer_(
        factory::make_default_load_balancer(std::forward<Args>(args)...)
      ) { };
  


  XCIntegratorImpl( const XCIntegratorImpl& )          = default;
  XCIntegratorImpl( XCIntegratorImpl&&      ) noexcept = default;


  virtual ~XCIntegratorImpl() noexcept = default;


  exc_vxc_type eval_exc_vxc( const MatrixType& P ) {
    return eval_exc_vxc_(P);
  }

};

}
}
