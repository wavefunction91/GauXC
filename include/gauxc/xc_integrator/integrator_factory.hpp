#pragma once

#include <gauxc/xc_integrator/integrator_defaults.hpp>

namespace GauXC  {
namespace detail {

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

      return make_default_host_integrator<MatrixType>(
#ifdef GAUXC_ENABLE_MPI
        comm,
#endif
        std::make_shared<functional_type>(func),
        std::make_shared<BasisSet<typename MatrixType::value_type>>(basis),
        lb
      );

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
#elif defined(GAUXC_ENABLE_SYCL)
      return make_default_sycl_integrator<MatrixType>(
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

}
}
