#pragma once

#include <gauxc/xc_integrator.hpp>

namespace GauXC  {
namespace detail {

/** Base class for XCIntegrator implementation */
template <typename MatrixType>
class XCIntegratorImpl {

public:

  using matrix_type    = MatrixType;
  using value_type     = typename matrix_type::value_type;
  using exc_vxc_type   = typename XCIntegrator<MatrixType>::exc_vxc_type;
  using exc_grad_type  = typename XCIntegrator<MatrixType>::exc_grad_type;
  using exx_type       = typename XCIntegrator<MatrixType>::exx_type;

protected:

  virtual exc_vxc_type  eval_exc_vxc_ ( const MatrixType& P ) = 0;
  virtual exc_grad_type eval_exc_grad_( const MatrixType& P ) = 0;
  virtual exx_type      eval_exx_     ( const MatrixType& P ) = 0;
  virtual const util::Timer& get_timings_() const = 0;
  virtual const LoadBalancer& get_load_balancer_() const = 0;
  
public:

  // Default all ctors as base is stateless

  XCIntegratorImpl()                                   = default;
  XCIntegratorImpl( const XCIntegratorImpl& )          = default;
  XCIntegratorImpl( XCIntegratorImpl&&      ) noexcept = default;
  virtual ~XCIntegratorImpl()                 noexcept = default;


  /** Integrate EXC / VXC (Mean field terms) for RKS
   * 
   *   TODO: add API for UKS/GKS
   *
   *  @param[in] P The alpha density matrix
   *  @returns EXC / VXC in a combined structure
   */
  exc_vxc_type eval_exc_vxc( const MatrixType& P ) {
    return eval_exc_vxc_(P);
  }

  /** Integrate EXC gradient for RKS
   * 
   *   TODO: add API for UKS/GKS
   *
   *  @param[in] P The alpha density matrix
   *  @returns EXC gradient
   */
  exc_grad_type eval_exc_grad( const MatrixType& P ) {
    return eval_exc_grad_(P);
  }

  /** Integrate Exact Exchange for RHF
   * 
   *   TODO: add API for UHF/GHF
   *
   *  @param[in] P The alpha density matrix
   *  @returns Excact Exchange Matrix
   */
  exx_type eval_exx( const MatrixType& P ) {
    return eval_exx_(P);
  }

  /** Get internal timers
   *
   *  @returns Timer instance for internal timings
   */
  const util::Timer& get_timings() const {
    return get_timings_();
  }


  const LoadBalancer& load_balancer() const {
    return get_load_balancer_();
  }
};

}
}
