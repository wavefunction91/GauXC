/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once

#include <gauxc/xc_integrator/pgas_dist/pgas_dist_xc_integrator_impl.hpp>

namespace GauXC {
namespace detail {

/// Base class for PGASDistributedXCIntegrator implentations on Host execution spaces
template <typename ValueType>
class PGASDistributedXCHostIntegrator : public PGASDistributedXCIntegratorImpl<ValueType> {

  using base_type  = PGASDistributedXCIntegratorImpl<ValueType>;

public:

  using value_type = typename base_type::value_type;
  using basis_type = typename base_type::basis_type;

  template <typename... Args>
  PGASDistributedXCHostIntegrator( Args&&... args) :
    base_type( std::forward<Args>(args)... ) { }

  virtual ~PGASDistributedXCHostIntegrator() noexcept;

};

extern template class PGASDistributedXCHostIntegrator<double>;



/// Factory to generate PGASDistributedXCHostIntegrator instances
template <typename ValueType>
struct PGASDistributedXCHostIntegratorFactory {

  using impl_type = PGASDistributedXCIntegratorImpl<ValueType>;
  using ptr_return_t = std::unique_ptr<impl_type>;

  /** Generate a PGASDistributedXCHostIntegrator instance
   *
   *  @param[in]  integration_kernel Name of integration scaffold to load ("Default", "Reference", etc)
   *  @param[in]  func               XC functional to integrate
   *  @param[in]  lb                 Pregenerated LoadBalancer instance
   *  @param[in]  lwd                Local Work Driver
   */
  static ptr_return_t make_integrator_impl( 
    std::string integrator_kernel,
    std::shared_ptr<functional_type>   func,
    std::shared_ptr<LoadBalancer>      lb,
    std::unique_ptr<LocalWorkDriver>&& lwd
    );

};


extern template struct PGASDistributedXCHostIntegratorFactory<double>;


}
}
