/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once

#include <memory>

#include <gauxc/types.hpp>
#include <gauxc/load_balancer.hpp>
#include <gauxc/xc_integrator_settings.hpp>

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

  using exc_vxc_type_rks  = std::tuple< value_type, matrix_type >;
  using exc_vxc_type_uks  = std::tuple< value_type, matrix_type, matrix_type >;  
  using exc_vxc_type_gks  = std::tuple< value_type, matrix_type, matrix_type, matrix_type, matrix_type >;
  using exc_grad_type = std::vector< value_type >;
  using exx_type      = matrix_type;

private:

  using pimpl_type    = detail::XCIntegratorImpl<MatrixType>;

  std::unique_ptr<pimpl_type> pimpl_;

public:

  XCIntegrator() = default;
  ~XCIntegrator() noexcept;

  XCIntegrator( std::unique_ptr<pimpl_type>&& pimpl );

  XCIntegrator( const XCIntegrator& ) = delete;
  XCIntegrator( XCIntegrator&& ) noexcept;

  value_type    integrate_den( const MatrixType& );

  value_type    eval_exc( const MatrixType&, const IntegratorSettingsXC& = IntegratorSettingsXC{} );
  value_type    eval_exc( const MatrixType&, const MatrixType&, const IntegratorSettingsXC& = IntegratorSettingsXC{} );
  value_type    eval_exc( const MatrixType&, const MatrixType&, const MatrixType&, const MatrixType&,  const IntegratorSettingsXC& = IntegratorSettingsXC{} );

  exc_vxc_type_rks  eval_exc_vxc ( const MatrixType&, 
                                   const IntegratorSettingsXC& = IntegratorSettingsXC{} );
  exc_vxc_type_uks  eval_exc_vxc ( const MatrixType&, const MatrixType&,
                                   const IntegratorSettingsXC& = IntegratorSettingsXC{} );
  exc_vxc_type_gks  eval_exc_vxc ( const MatrixType&, const MatrixType&, const MatrixType&, const MatrixType&,
                                   const IntegratorSettingsXC& = IntegratorSettingsXC{});

  exc_grad_type eval_exc_grad( const MatrixType& );

  exx_type      eval_exx     ( const MatrixType&, 
                               const IntegratorSettingsEXX& = IntegratorSettingsEXX{} );


  const util::Timer& get_timings() const;
  const LoadBalancer& load_balancer() const;
  LoadBalancer& load_balancer();
};


}
