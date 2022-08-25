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

  using exc_vxc_type  = std::tuple< value_type, matrix_type >;
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
  exc_vxc_type  eval_exc_vxc ( const MatrixType& );
  exc_grad_type eval_exc_grad( const MatrixType& );
  exx_type      eval_exx     ( const MatrixType&, 
                               const IntegratorSettingsEXX& = IntegratorSettingsEXX{} );


  const util::Timer& get_timings() const;
  const LoadBalancer& load_balancer() const;
};


}
