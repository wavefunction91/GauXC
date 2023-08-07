/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once

#include <gauxc/xc_integrator/replicated_xc_integrator.hpp>
#include <gauxc/xc_integrator/local_work_driver.hpp>
#include <gauxc/reduction_driver.hpp>
#include <gauxc/types.hpp>
#include <gauxc/basisset.hpp>

namespace GauXC  {
namespace detail {


/// Base class for ReplicatedXCIntegrator implementations
template <typename ValueType>
class ReplicatedXCIntegratorImpl {

public:

  using value_type = ValueType;
  using basis_type = BasisSet< value_type >;

protected:

  std::shared_ptr< functional_type > func_;               ///< XC functional
  std::shared_ptr< LoadBalancer >    load_balancer_;      ///< Load Balancer
  std::unique_ptr< LocalWorkDriver > local_work_driver_;  ///< Local Work Driver
  std::shared_ptr< ReductionDriver > reduction_driver_;   ///< Reduction Driver

  util::Timer timer_;


  virtual void integrate_den_( int64_t m, int64_t n, const value_type* P,
                               int64_t ldp, value_type* N_EL ) = 0;
  virtual void eval_exc_vxc_( int64_t m, int64_t n, const value_type* P,
                              int64_t ldp, value_type* VXC, int64_t ldvxc,
                              value_type* EXC ) = 0;
  virtual void eval_exc_vxc_( int64_t m, int64_t n, const value_type* Pscalar,
                              int64_t ldpscalar,
                              const value_type* Pz,
                              int64_t ldpz,
                              value_type* VXCscalar, int64_t ldvxcscalar,
                              value_type* VXCz, int64_t ldvxcz,
                              value_type* EXC ) = 0;
  virtual void eval_exc_grad_( int64_t m, int64_t n, const value_type* P,
                               int64_t ldp, value_type* EXC_GRAD ) = 0;
  virtual void eval_exx_( int64_t m, int64_t n, const value_type* P,
                          int64_t ldp, value_type* K, int64_t ldk,
                          const IntegratorSettingsEXX& settings ) = 0;

public:

  ReplicatedXCIntegratorImpl( std::shared_ptr< functional_type >   func,
                              std::shared_ptr< LoadBalancer >      lb, 
                              std::unique_ptr< LocalWorkDriver >&& lwd,
                              std::shared_ptr< ReductionDriver>    rd
                              );

  virtual ~ReplicatedXCIntegratorImpl() noexcept;

  void integrate_den( int64_t m, int64_t n, const value_type* P,
                      int64_t ldp, value_type* N_EL );

  void eval_exc_vxc( int64_t m, int64_t n, const value_type* P,
                     int64_t ldp, value_type* VXC, int64_t ldvxc,
                     value_type* EXC ); 

  void eval_exc_vxc( int64_t m, int64_t n, const value_type* Pscalar,
                     int64_t ldpscalar,
                     const value_type* Pz,
                     int64_t ldpz,
                     value_type* VXCscalar, int64_t ldvxcscalar,
                     value_type* VXCz, int64_t ldvxcz,
                     value_type* EXC );

  void eval_exc_grad( int64_t m, int64_t n, const value_type* P,
                      int64_t ldp, value_type* EXC_GRAD );

  void eval_exx( int64_t m, int64_t n, const value_type* P,
                 int64_t ldp, value_type* K, int64_t ldk,
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


extern template class ReplicatedXCIntegratorImpl<double>;

}
}
