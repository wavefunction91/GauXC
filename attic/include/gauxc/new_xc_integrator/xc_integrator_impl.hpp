#pragma once

#include <gauxc/xc_integrator.hpp>

namespace GauXC  {
namespace detail {

template <typename MatrixType>
class XCIntegratorImpl {

public:

  using matrix_type   = MatrixType;
  using value_type    = typename matrix_type::value_type;
  using exc_vxc_type  = typename XCIntegrator<MatrixType>::exc_vxc_type;

protected:

  virtual exc_vxc_type eval_exc_vxc_( const MatrixType& ) = 0;
  virtual const util::Timer& get_timings_() const = 0;
  
public:

  XCIntegratorImpl()                                   = default;
  XCIntegratorImpl( const XCIntegratorImpl& )          = default;
  XCIntegratorImpl( XCIntegratorImpl&&      ) noexcept = default;
  virtual ~XCIntegratorImpl()                 noexcept = default;


  exc_vxc_type eval_exc_vxc( const MatrixType& P ) {
    return eval_exc_vxc_(P);
  }

  const util::Timer& get_timings() const {
    return get_timings_();
  }

};

}
}
