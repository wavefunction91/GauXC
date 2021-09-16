#include "hip_backend.hpp"

namespace GauXC {

HIPBackend::HIPBackend() {

  // Create HIP Stream and CUBLAS Handles and make them talk to eachother
  master_stream = std::make_shared< util::hip_stream >();
  master_handle = std::make_shared< util::hipblas_handle >();

  hipblasSetStream( *master_handle, *master_stream );

}

HIPBackend::~HIPBackend() noexcept = default;

HIPBackend::device_buffer_t HIPBackend::allocate_device_buffer(int64_t sz) {
  void* ptr;
  auto stat = hipMalloc(&ptr, sz);
  GAUXC_HIP_ERROR( "HIP Malloc Failed", stat );
  return device_buffer_t{ptr,sz};
}

size_t HIPBackend::get_available_mem() {
  size_t hip_avail, hip_total;
  auto stat = hipMemGetInfo( &hip_avail, &hip_total );
  GAUXC_HIP_ERROR( "MemInfo Failed", stat );
  return hip_avail;
}

void HIPBackend::free_device_buffer( void* ptr ) {
  auto stat = hipFree(ptr);
  GAUXC_HIP_ERROR( "Free Failed", stat );
}

void HIPBackend::master_queue_synchronize() {
  auto stat = hipStreamSynchronize( *master_stream );
  GAUXC_HIP_ERROR( "StreamSynchronized Failed", stat );
}

type_erased_queue HIPBackend::queue() {
  return type_erased_queue(master_stream);
}

void HIPBackend::create_blas_queue_pool(int32_t ns) {
  blas_streams.resize(ns);
  blas_handles.resize(ns);
  for( auto i = 0; i < ns; ++i )
    hipblasSetStream( blas_handles[i], blas_streams[i] );
}

void HIPBackend::copy_async_( size_t sz, const void* src, void* dest,
  std::string msg ) {
  auto stat = hipMemcpyAsync( dest, src, sz, hipMemcpyDefault, *master_stream );
  GAUXC_HIP_ERROR( "HIP Memcpy Async Failed ["+msg+"]", stat );
}

void HIPBackend::set_zero_(size_t sz, void* data, std::string msg ) {
  auto stat = hipMemset( data, 0, sz );
  GAUXC_HIP_ERROR( "HIP Memset Failed ["+msg+"]", stat );
}

void HIPBackend::copy_async_2d_( size_t M, size_t N, const void* A, size_t LDA,
  void* B, size_t LDB, std::string msg ) {
  auto stat = hipMemcpy2DAsync( B, LDB, A, LDA, M, N, hipMemcpyDefault,
    *master_stream );
  GAUXC_HIP_ERROR( "HIP 2D Memcpy Async Failed ["+msg+"]", stat );
}

}
