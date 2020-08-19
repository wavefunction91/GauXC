#include <gauxc/xc_integrator/xc_cuda_data.hpp>
#include "buffer_adaptor.hpp"

namespace GauXC {

template <typename F>
XCCudaData<F>::XCCudaData( size_t _natoms,
                           size_t _n_deriv, 
                           size_t _nbf,
                           size_t _nshells,
                           bool _denpack_host,
                           bool _vxcinc_host  ):
  nshells(_nshells), 
  nbf(_nbf), 
  n_deriv(_n_deriv), 
  natoms(_natoms),
  denpack_host(_denpack_host), 
  vxcinc_host(_vxcinc_host)  {


  // TODO: Expose this
  double fill_fraction = 0.9;

  cudaError_t stat;

  // Get Total Available Memory
  size_t cuda_avail, cuda_total;
  stat = cudaMemGetInfo( &cuda_avail, &cuda_total );
  GAUXC_CUDA_ERROR( "MemInfo Failed", stat );

  // Allocate up to fill_fraction
  size_t fill_sz = fill_fraction * cuda_total;
  stat = cudaMalloc( &device_ptr, fill_sz );
  GAUXC_CUDA_ERROR( "CUDA Malloc Failed", stat );

  // Allocate static memory with proper alignment
  buffer_adaptor mem( device_ptr, fill_sz );


  shells_device = mem.aligned_alloc<Shell<F>>( nshells );
  exc_device    = mem.aligned_alloc<F>( 1 );
  nel_device    = mem.aligned_alloc<F>( 1 );
  rab_device    = mem.aligned_alloc<F>( natoms * natoms );
  coords_device = mem.aligned_alloc<F>( 3 * natoms );

  if( not vxcinc_host )
    vxc_device = mem.aligned_alloc<F>( nbf * nbf );
  if( not denpack_host )
    dmat_device = mem.aligned_alloc<F>( nbf * nbf );

  // Get current stack location
  dynmem_ptr = mem.stack();
  dynmem_sz  = mem.nleft(); 

  // Create CUDA Stream and CUBLAS Handles and make them talk to eachother
  master_stream = std::make_unique< util::cuda_stream >();
  master_handle = std::make_unique< util::cublas_handle >();

  cublasSetStream( *master_handle, *master_stream );

  // Create MAGMA Queue from CUDA Stream and CUBLAS Handle
  master_magma_queue = 
    std::make_unique< util::magma_queue >( 0, *master_stream, *master_handle );

}



template <typename F>
XCCudaData<F>::~XCCudaData() noexcept {
  if( device_ptr ) util::cuda_free( device_ptr );
} 




// Explicit Instantiations
template struct XCCudaData<double>;

}
