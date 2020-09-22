#pragma once

#include <gauxc/xc_integrator/default_xc_host_integrator.hpp>
#include <gauxc/xc_integrator/default_xc_cuda_integrator.hpp>
#include <gauxc/xc_integrator/default_xc_hip_integrator.hpp>
#include <gauxc/xc_integrator/default_xc_sycl_integrator.hpp>

namespace GauXC  {
namespace detail {

template <typename MatrixType, typename... Args>
std::unique_ptr<XCIntegratorImpl<MatrixType>> make_default_host_integrator( Args&&... args ) {
  return std::make_unique<DefaultXCHostIntegrator<MatrixType>>( std::forward<Args>(args)... );
}

#ifdef GAUXC_ENABLE_CUDA
template <typename MatrixType, typename... Args>
std::unique_ptr<XCIntegratorImpl<MatrixType>> make_default_cuda_integrator( Args&&... args ) {
  return std::make_unique<DefaultXCCudaIntegrator<MatrixType>>( std::forward<Args>(args)... );
}
#endif

#ifdef GAUXC_ENABLE_SYCL
template <typename MatrixType, typename... Args>
std::unique_ptr<XCIntegratorImpl<MatrixType>> make_default_sycl_integrator( Args&&... args ) {
  return std::make_unique<DefaultXCSyclIntegrator<MatrixType>>( std::forward<Args>(args)... );
}
#endif

#ifdef GAUXC_ENABLE_HIP
template <typename MatrixType, typename... Args>
std::unique_ptr<XCIntegratorImpl<MatrixType>> make_default_hip_integrator( Args&&... args ) {
  return std::make_unique<DefaultXCHipIntegrator<MatrixType>>( std::forward<Args>(args)... );
}
#endif

}
}
