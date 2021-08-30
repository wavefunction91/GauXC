#pragma once

#include <gauxc/new_xc_integrator/xc_integrator_impl.hpp>

namespace GauXC  {
namespace detail {

template <typename ValueType>
class ReplicatedXCIntegratorImpl;

template <typename MatrixType>
class ReplicatedXCIntegrator : public XCIntegratorImpl<MatrixType> {

public:

  using matrix_type   = typename XCIntegratorImpl<MatrixType>::matrix_type;
  using value_type    = typename XCIntegratorImpl<MatrixType>::value_type;
  using exc_vxc_type  = typename XCIntegratorImpl<MatrixType>::exc_vxc_type;

private:

  using pimpl_type = ReplicatedXCIntegratorImpl<value_type>;
  std::unique_ptr< pimpl_type > pimpl_;

  exc_vxc_type eval_exc_vxc_( const MatrixType& ) override;
  const util::Timer& get_timings_() const override;

public:

  ReplicatedXCIntegrator();
  ReplicatedXCIntegrator( std::unique_ptr<pimpl_type>&& );

  ~ReplicatedXCIntegrator() noexcept;

  ReplicatedXCIntegrator( const ReplicatedXCIntegrator& ) = delete;
  ReplicatedXCIntegrator( ReplicatedXCIntegrator&& ) noexcept;

};


}
}
