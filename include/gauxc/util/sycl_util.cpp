#pragma once

#include <CL/sycl.hpp>
#include <string>
#include <gauxc/gauxc_config.hpp>

#ifdef GAUXC_ENABLE_SYCL

namespace GauXC {
namespace util  {

template <typename T>
inline T* sycl_malloc( size_t n, cl::sycl::queue& q ) {
    T* ptr = cl::sycl::malloc_device<T>(n * sizeof(T), q);
    return ptr;
}

template <typename T>
inline void sycl_free( T*& ptr, cl::sycl::queue& q ) {
  cl::sycl::free((void *)ptr, q);
  ptr = nullptr;
}

template <typename T>
inline void sycl_copy( size_t len, T* dest, const T* src, cl::sycl::queue& q, std::string m = "") {
  auto stat = q.memcpy(dest, src, len * sizeof(T));
  q.wait();
}
template <typename T>
inline void sycl_copy_async(size_t len, T *dest, const T *src, cl::sycl::queue& q, std::string m = "") {
  auto stat = q.memcpy(dest, src, len * sizeof(T));
}


template <typename T>
inline void sycl_set_zero( size_t len, T* ptr, cl::sycl::queue& q, std::string m = "") {
  auto stat = q.memset(ptr, 0, len * sizeof(T));
  q.wait()
}
template <typename T>
inline void sycl_set_zero_async(size_t len, T *ptr, cl::sycl::queue& q, std::string m = "") {
  auto stat = q.memset(ptr, 0, len * sizeof(T));
}


inline void sycl_device_sync(cl::sycl::queue& q) {
  q.wait_and_throw();
}

}
}

#endif
