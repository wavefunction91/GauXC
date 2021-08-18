#include "xc_cuda_aos_data_base.hpp"

namespace GauXC {

XCCudaAoSDataBase::XCCudaAoSDataBase() {

  // Create CUDA Stream and CUBLAS Handles and make them talk to eachother
  master_stream = std::make_unique< util::cuda_stream >();
  master_handle = std::make_unique< util::cublas_handle >();

  cublasSetStream( *master_handle, *master_stream );

}

XCCudaAoSDataBase::~XCCudaAoSDataBase() noexcept = default;


void XCCudaAoSDataBase::allocate_device_buffer(int64_t sz) {

  if( sz > 0 ) throw std::runtime_error("Fixed Size Device Allocations NYI");

  // TODO: Expose this
  double fill_fraction = 0.9;

  cudaError_t stat;

  // Get Total Available Memory
  size_t cuda_avail, cuda_total;
  stat = cudaMemGetInfo( &cuda_avail, &cuda_total );
  GAUXC_CUDA_ERROR( "MemInfo Failed", stat );

  // Allocate up to fill_fraction
  devmem_sz = fill_fraction * cuda_avail;
  stat = cudaMalloc( &device_ptr, devmem_sz );
  GAUXC_CUDA_ERROR( "CUDA Malloc Failed", stat );

}

void XCCudaAoSDataBase::master_queue_synchronize() {
  cudaStreamSynchronize( *master_stream );
}

void XCCudaAoSDataBase::copy_async_( size_t sz, const void* src, void* dest,
  std::string msg ) {
  util::cuda_copy_async( sz, (char*)dest, (const char*)src, *master_stream, msg );
}

void XCCudaAoSDataBase::set_zero_(size_t sz, void* data, std::string msg ) {
  util::cuda_set_zero(sz, (char*)data, msg );
}

}
