#pragma once

#include <gauxc/xc_integrator/integrator_defaults.hpp>

namespace GauXC  {
namespace detail {

template <typename MatrixType>
std::unique_ptr<XCIntegratorImpl<MatrixType>> 
  integrator_factory( ExecutionSpace ex, MPI_Comm comm, const functional_type& func, 
                      const BasisSet<typename MatrixType::value_type>& basis, 
                      std::shared_ptr<LoadBalancer> lb) {

    if( ex == ExecutionSpace::Host ) {

#ifdef GAUXC_ENABLE_HOST

      return make_default_host_integrator<MatrixType>(
        comm,
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
        comm,
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

}
}
