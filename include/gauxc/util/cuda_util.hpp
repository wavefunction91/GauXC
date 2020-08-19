#pragma once
#include <gauxc/gauxc_config.hpp>
#include <gauxc/exceptions/cuda_exception.hpp>

#ifdef GAUXC_ENABLE_CUDA

namespace GauXC {
namespace util  {

struct cuda_stream {

  cudaStream_t stream;
  inline cuda_stream() {
    auto stat = cudaStreamCreate( &stream );
    GAUXC_CUDA_ERROR("CUDA Stream Create Failed", stat);
  }

  inline ~cuda_stream() noexcept {
    if( stream != 0 ) cudaStreamDestroy( stream );
  }

  cuda_stream( const cuda_stream& ) = delete;
  inline cuda_stream( cuda_stream&& other ) noexcept {
    stream = other.stream;
    other.stream = 0;
  };

  inline operator cudaStream_t() const { return stream; }

};








template <typename T>
inline T* cuda_malloc( size_t n ) {

  T* ptr;
  auto stat = cudaMalloc( (void**)&ptr, n * sizeof(T) );
  GAUXC_CUDA_ERROR( "CUDA Malloc Failed", stat );

  return ptr;
}






template <typename T>
inline void cuda_free( T*& ptr ) {
  auto stat = cudaFree( (void*)ptr );
  GAUXC_CUDA_ERROR( "CUDA Free Failed", stat );
  ptr = nullptr;
}

template <typename T, typename... Args>
inline void cuda_free( T*& ptr, Args&&... args ) {
  cuda_free(ptr);
  cuda_free(std::forward<Args>(args)...);
}





template <typename T>
inline void cuda_copy( size_t len, T* dest, const T* src ) {
  auto stat = cudaMemcpy( dest, src, len * sizeof(T), cudaMemcpyDefault );
  GAUXC_CUDA_ERROR( "CUDA Memcpy Failed", stat );
}

template <typename T>
inline void cuda_copy_async( size_t len, T* dest, const T* src, cudaStream_t s ) {
  auto stat = cudaMemcpyAsync( dest, src, len * sizeof(T), cudaMemcpyDefault, s );
  GAUXC_CUDA_ERROR( "CUDA Memcpy Async Failed", stat );
}




inline void cuda_device_sync() {
  auto stat = cudaDeviceSynchronize();
  GAUXC_CUDA_ERROR( "CUDA Device Sync Failed", stat );
}

}
}

#endif
