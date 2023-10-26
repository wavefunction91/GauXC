/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once

#include <gauxc/xc_integrator/pgas_dist_xc_integrator.hpp>
#include <gauxc/xc_integrator/local_work_driver.hpp>
#include <gauxc/dist_array.hpp>
#include <gauxc/reduction_driver.hpp>
#include <gauxc/types.hpp>
#include <gauxc/basisset.hpp>

namespace GauXC  {
namespace detail {


/// Base class for PGASDistributedXCIntegrator implementations
template <typename ValueType>
class PGASDistributedXCIntegratorImpl {

public:

  using value_type = ValueType;
  using matrix_type = Darray; // TODO: template on value_type
  using basis_type = BasisSet< value_type >;

protected:

  std::shared_ptr< functional_type > func_;               ///< XC functional
  std::shared_ptr< LoadBalancer >    load_balancer_;      ///< Load Balancer
  std::unique_ptr< LocalWorkDriver > local_work_driver_;  ///< Local Work Driver

  util::Timer timer_;


  virtual void integrate_den_( const matrix_type& P, value_type* N_EL ) = 0;
  virtual void eval_exc_vxc_( const matrix_type& P, matrix_type& VXC, value_type* EXC ) = 0;
  virtual void eval_exc_vxc_( const matrix_type& Ps, const matrix_type& Pz, 
                              matrix_type& VXCs, matrix_type& VXCz, value_type* EXC ) = 0;
  virtual void eval_exc_grad_( const matrix_type& P, value_type* EXC_GRAD ) = 0;
  virtual void eval_exx_( const matrix_type& P, matrix_type& K,
                          const IntegratorSettingsEXX& settings ) = 0;

public:

  PGASDistributedXCIntegratorImpl( std::shared_ptr< functional_type >   func,
                              std::shared_ptr< LoadBalancer >      lb, 
                              std::unique_ptr< LocalWorkDriver >&& lwd
                              );

  virtual ~PGASDistributedXCIntegratorImpl() noexcept;

  void integrate_den( const matrix_type& P, value_type* N_EL );
  void eval_exc_vxc( const matrix_type& P, matrix_type& VXC, value_type* EXC );
  void eval_exc_vxc( const matrix_type& Ps, const matrix_type& Pz, 
                     matrix_type& VXCs, matrix_type& VXCz, value_type* EXC );
  void eval_exc_grad( const matrix_type& P, value_type* EXC_GRAD );
  void eval_exx( const matrix_type& P, matrix_type& K,
                 const IntegratorSettingsEXX& settings );

  inline const util::Timer& get_timings() const { return timer_; }

  inline std::unique_ptr< LocalWorkDriver > release_local_work_driver() {
    return std::move( local_work_driver_ );
  }

  inline const auto& load_balancer() const { return *load_balancer_; }
  inline auto& load_balancer() { return *load_balancer_; }
  inline const auto& get_load_balancer() const { return load_balancer(); }
  inline auto& get_load_balancer() { return load_balancer(); }
};


extern template class PGASDistributedXCIntegratorImpl<double>;

}
}
