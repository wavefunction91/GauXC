#pragma once
#include <gauxc/gauxc_config.hpp>
#include <gauxc/exceptions/hip_exception.hpp>

#ifdef GAUXC_ENABLE_HIP

namespace GauXC {
namespace util  {

struct hip_stream;
struct hip_event;

struct hip_stream {

  hipStream_t stream;
  inline hip_stream() {
    auto stat = hipStreamCreate( &stream );
    GAUXC_HIP_ERROR("HIP Stream Create Failed", stat);
  }

  inline ~hip_stream() noexcept {
    if( stream != 0 ) hipStreamDestroy( stream );
  }

  hip_stream( const hip_stream& ) = delete;
  inline hip_stream( hip_stream&& other ) noexcept {
    stream = other.stream;
    other.stream = 0;
  };

  inline operator hipStream_t() const { return stream; }

  inline void wait( hipEvent_t event ) {
    auto stat = hipStreamWaitEvent( stream, event, 0 );
    GAUXC_HIP_ERROR("STREAM WAIT FAILED", stat );
  }
};


struct hip_event {

  hipEvent_t event;
  inline hip_event() {
    auto stat = hipEventCreate( &event );
    GAUXC_HIP_ERROR("HIP Event Create Failed", stat);
  }

  inline ~hip_event() noexcept {
    if( event != 0 ) hipEventDestroy( event );
  }

  hip_event( const hip_event& ) = delete;
  inline hip_event( hip_event&& other ) noexcept {
    event = other.event;
    other.event = 0;
  };

  inline operator hipEvent_t() const { return event; }

  inline void record( hipStream_t stream ) {
    auto stat = hipEventRecord( event, stream );
    GAUXC_HIP_ERROR("Event Record Failed", stat );
  }

};





template <typename T>
inline T* hip_malloc( size_t n ) {

  T* ptr;
  auto stat = hipMalloc( (void**)&ptr, n * sizeof(T) );
  GAUXC_HIP_ERROR( "HIP Malloc Failed", stat );

  return ptr;
}






template <typename T>
inline void hip_free( T*& ptr ) {
  auto stat = hipFree( (void*)ptr );
  GAUXC_HIP_ERROR( "HIP Free Failed", stat );
  ptr = nullptr;
}

template <typename T, typename... Args>
inline void hip_free( T*& ptr, Args&&... args ) {
  hip_free(ptr);
  hip_free(std::forward<Args>(args)...);
}





template <typename T>
inline void hip_copy( size_t len, T* dest, const T* src, std::string m = "") {
  auto stat = hipMemcpy( dest, src, len * sizeof(T), hipMemcpyDefault );
  GAUXC_HIP_ERROR( "HIP Memcpy Failed ["+m+"]", stat );
}

template <typename T>
inline void hip_copy_async( size_t len, T* dest, const T* src, hipStream_t s,
                             std::string m = "" ) {
  auto stat = hipMemcpyAsync( dest, src, len * sizeof(T), hipMemcpyDefault, s );
  GAUXC_HIP_ERROR( "HIP Memcpy Async Failed ["+m+"]", stat );
}


template <typename T>
inline void hip_set_zero( size_t len, T* ptr, std::string m = "" ) {
  auto stat = hipMemset( ptr, 0, len * sizeof(T) );
  GAUXC_HIP_ERROR( "HIP Memset Failed ["+m+"]", stat );
}

template <typename T>
inline void hip_set_zero_async( size_t len, T* ptr, hipStream_t stream, 
                                 std::string m = "" ) {
  auto stat = hipMemsetAsync( ptr, 0, len * sizeof(T), stream );
  GAUXC_HIP_ERROR( "HIP Memset Async Failed ["+m+"]", stat );
}



inline void hip_device_sync() {
  auto stat = hipDeviceSynchronize();
  GAUXC_HIP_ERROR( "HIP Device Sync Failed", stat );
}



template <typename T>
inline int hip_kernel_max_threads_per_block( T* func ) {
  hipFuncAttributes attr;
  auto stat = hipFuncGetAttributes(&attr, (void*)func);

  GAUXC_HIP_ERROR( "GetAttr Failed", stat ); 
  return attr.maxThreadsPerBlock;
}

}
}

#endif
