#include "cuda_backend.hpp"

namespace GauXC {

CUDABackend::CUDABackend() {

  // Create CUDA Stream and CUBLAS Handles and make them talk to eachother
  master_stream = std::make_unique< util::cuda_stream >();
  master_handle = std::make_unique< util::cublas_handle >();

  cublasSetStream( *master_handle, *master_stream );

}

CUDABackend::~CUDABackend() noexcept = default;

CUDABackend::device_buffer_t CUDABackend::allocate_device_buffer(int64_t sz) {
  void* ptr;
  auto stat = cudaMalloc(&ptr, sz);
  GAUXC_CUDA_ERROR( "CUDA Malloc Failed", stat );
  return device_buffer_t{ptr,sz};
}

size_t CUDABackend::get_available_mem() {
  size_t cuda_avail, cuda_total;
  auto stat = cudaMemGetInfo( &cuda_avail, &cuda_total );
  GAUXC_CUDA_ERROR( "MemInfo Failed", stat );
  return cuda_avail;
}

void CUDABackend::free_device_buffer( void* ptr ) {
  cudaFree(ptr);
}

void CUDABackend::master_queue_synchronize() {
  cudaStreamSynchronize( *master_stream );
}

void CUDABackend::create_blas_queue_pool(int32_t ns) {
  blas_streams.resize(ns);
  blas_handles.resize(ns);
  for( auto i = 0; i < ns; ++i )
    cublasSetStream( blas_handles[i], blas_streams[i] );
}

void CUDABackend::copy_async_( size_t sz, const void* src, void* dest,
  std::string msg ) {
  auto stat = cudaMemcpyAsync( dest, src, sz, cudaMemcpyDefault, *master_stream );
  GAUXC_CUDA_ERROR( "CUDA Memcpy Async Failed ["+msg+"]", stat );
}

void CUDABackend::set_zero_(size_t sz, void* data, std::string msg ) {
  auto stat = cudaMemset( data, 0, sz );
  GAUXC_CUDA_ERROR( "CUDA Memset Failed ["+msg+"]", stat );
}

void CUDABackend::copy_async_2d_( size_t M, size_t N, const void* A, size_t LDA,
  void* B, size_t LDB, std::string msg ) {
  auto stat = cudaMemcpy2DAsync( B, LDB, A, LDA, M, N, cudaMemcpyDefault,
    *master_stream );
  GAUXC_CUDA_ERROR( "CUDA 2D Memcpy Async Failed ["+msg+"]", stat );
}

}
