#pragma once

#include <gauxc/xc_integrator/default_xc_host_integrator.hpp>

namespace GauXC  {
namespace detail {

template <typename MatrixType, typename... Args>
std::unique_ptr<XCIntegratorImpl<MatrixType>> make_default_host_integrator( Args&&... args ) {
  return std::make_unique<DefaultXCHostIntegrator<MatrixType>>( std::forward<Args>(args)... );
}

}
}
