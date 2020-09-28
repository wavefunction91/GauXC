#pragma once
#include <gauxc/gauxc_config.hpp>
#include <gauxc/exceptions/magma_exception.hpp>

#ifdef GAUXC_ENABLE_CUDA
#ifdef GAUXC_ENABLE_MAGMA

namespace GauXC {
namespace util  {

struct magma_queue {

  magma_queue_t queue;

  

  inline magma_queue(magma_int_t dev) {
    magma_queue_create( dev, &queue );
  }

  inline magma_queue() : magma_queue(0) { }

  inline magma_queue( magma_int_t dev, cudaStream_t stream, cublasHandle_t handle ) {
    magma_queue_create_from_cuda( dev, stream, handle, NULL, &queue );
  }

  inline magma_queue( cudaStream_t stream, cublasHandle_t handle ) :
    magma_queue( 0, stream, handle ) { }

  inline ~magma_queue() noexcept {
    if( queue != 0 ) magma_queue_destroy( queue );
  }

  magma_queue( const magma_queue& ) = delete;
  inline magma_queue( magma_queue&& other ) noexcept {
    queue = other.queue;
    other.queue = 0;
  };

  inline operator magma_queue_t() const { return queue; }

};

}
}

#endif
#endif
