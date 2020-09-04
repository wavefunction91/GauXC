#pragma once
#include <gauxc/gauxc_config.hpp>
#include <gauxc/exceptions/cublas_exception.hpp>

#ifdef GAUXC_ENABLE_CUDA

namespace GauXC {
namespace util  {

struct cublas_handle {

  cublasHandle_t handle;
  inline cublas_handle() {
    auto stat = cublasCreate( &handle );
    GAUXC_CUBLAS_ERROR("CUBLAS Handle Create Failed", stat);
  }

  inline ~cublas_handle() noexcept {
    if( handle != 0 ) cublasDestroy( handle );
  }

  cublas_handle( const cublas_handle& ) = delete;
  inline cublas_handle( cublas_handle&& other ) noexcept {
    handle = other.handle;
    other.handle = 0;
  };

  inline operator cublasHandle_t() const { return handle; }

};


inline static cudaStream_t get_stream( cublasHandle_t handle ) {
  cudaStream_t stream;
  auto stat = cublasGetStream(handle, &stream );  
  GAUXC_CUBLAS_ERROR("CUBLAS GET STREAM FAILED", stat );
  return stream;
}

}
}

#endif
