#pragma once

#include <gauxc/xc_integrator/reference_xc_host_integrator.hpp>
#include <gauxc/xc_integrator/incore_xc_cuda_integrator.hpp>
#include <gauxc/xc_integrator/shellbatched_xc_cuda_integrator.hpp>

namespace GauXC  {
namespace detail {

template <typename MatrixType>
using DefaultXCHostIntegrator = ReferenceXCHostIntegrator<MatrixType>;

template <typename MatrixType>
//using DefaultXCCudaIntegrator = IncoreXCCudaIntegrator<MatrixType>;
using DefaultXCCudaIntegrator = ShellBatchedXCCudaIntegrator<MatrixType>;


#ifdef GAUXC_ENABLE_HOST
template <typename MatrixType, typename... Args>
std::unique_ptr<XCIntegratorImpl<MatrixType>> 
  make_default_host_integrator( Args&&... args ) {
    return std::make_unique<DefaultXCHostIntegrator<MatrixType>>( 
      std::forward<Args>(args)... 
    );
}
#endif

#ifdef GAUXC_ENABLE_CUDA
template <typename MatrixType, typename... Args>
std::unique_ptr<XCIntegratorImpl<MatrixType>> 
  make_default_cuda_integrator( Args&&... args ) {
  return std::make_unique<DefaultXCCudaIntegrator<MatrixType>>( 
    std::forward<Args>(args)...
  );
}
#endif

}
}
