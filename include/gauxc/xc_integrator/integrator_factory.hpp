#pragma once

#include <gauxc/xc_integrator/integrator_defaults.hpp>
#include <gauxc/util/forward_as_shared_ptr.hpp>

namespace GauXC  {
namespace detail {


template <typename MatrixType, typename... Args>
std::unique_ptr<XCIntegratorImpl<MatrixType>> 
  default_integrator_factory( ExecutionSpace ex, Args&&... args ) {
  
  if( ex == ExecutionSpace::Host ) {

#ifdef GAUXC_ENABLE_HOST
    return make_default_host_integrator<MatrixType>( 
      forward_as_shared_ptr(args)... 
    );
#else
    throw std::runtime_error("GAUXC_ENABLE_HOST is FALSE");
    return nullptr;
#endif

  } else {

#ifdef GAUXC_ENABLE_CUDA
    return make_default_cuda_integrator<MatrixType>( forward_as_shared_ptr(args)... );
#else
    throw std::runtime_error("GAUXC_ENABLE_DEVICE is FALSE");
    return nullptr;
#endif

  }
}

}


template <typename MatrixType, typename... Args>
XCIntegrator<MatrixType>
  make_default_integrator( ExecutionSpace ex, Args&&... args ) {

  return XCIntegrator<MatrixType>(
    detail::default_integrator_factory<MatrixType>( ex, 
      std::forward<Args>(args)... 
    )
  );

}


}
