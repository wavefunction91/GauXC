#pragma once

#include <gauxc/xc_integrator/integrator_defaults.hpp>
#include <gauxc/util/forward_as_shared_ptr.hpp>

namespace GauXC  {
namespace detail {

#if 0
#ifdef GAUXC_ENABLE_MPI
template <typename MatrixType>
std::unique_ptr<XCIntegratorImpl<MatrixType>> 
  integrator_factory( ExecutionSpace ex, MPI_Comm comm, const functional_type& func, 
                      const BasisSet<typename MatrixType::value_type>& basis, 
                      std::shared_ptr<LoadBalancer> lb) 
#else
template <typename MatrixType>
std::unique_ptr<XCIntegratorImpl<MatrixType>> 
  integrator_factory( ExecutionSpace ex, const functional_type& func, 
                      const BasisSet<typename MatrixType::value_type>& basis, 
                      std::shared_ptr<LoadBalancer> lb) 
#endif
                      
{

    if( ex == ExecutionSpace::Host ) {

#ifdef GAUXC_ENABLE_HOST

      return make_default_host_integrator<MatrixType>(
#ifdef GAUXC_ENABLE_MPI
        comm,
#endif
        std::make_shared<functional_type>(func),
        std::make_shared<BasisSet<typename MatrixType::value_type>>(basis),
        lb
      );
#else

      throw std::runtime_error("GAUXC_ENABLE_HOST is FALSE");
      return nullptr;

#endif

    } else {

#ifdef GAUXC_ENABLE_CUDA

      return make_default_cuda_integrator<MatrixType>(
#ifdef GAUXC_ENABLE_MPI
        comm,
#endif
        std::make_shared<functional_type>(func),
        std::make_shared<BasisSet<typename MatrixType::value_type>>(basis),
        lb
      );

#else
      throw std::runtime_error("GAUXC_ENABLE_DEVICE is FALSE");
      return nullptr;
#endif

    }

  }

#else

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

#endif
}


template <typename MatrixType, typename... Args>
XCIntegrator<MatrixType>
  make_default_integrator( ExecutionSpace ex, Args&&... args ) {

  return XCIntegrator<MatrixType>(
    detail::default_integrator_factory<MatrixType>( ex, std::forward<Args>(args)... )
  );

}


}
