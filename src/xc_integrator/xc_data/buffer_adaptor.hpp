#pragma once
#include <memory>

//#define csl __PRETTY_FUNCTION__
#define csl std::string(__FILE__) + ": " + std::to_string(__LINE__)

namespace GauXC {

class buffer_adaptor {

  size_t nalloc_;
  size_t nleft_;
  void*  top_;
  void*  stack_;

public:

  buffer_adaptor() = delete;

  inline buffer_adaptor( void* ptr, size_t len ) :
    nalloc_(len), 
    nleft_(len), 
    top_(ptr), 
    stack_(ptr) { }

  template <typename T>
  T* aligned_alloc( size_t len, 
                    size_t align = alignof(T),
                    std::string msg = "" ) {

    char* old_stack = (char*)stack_;
    if( std::align( align, 
                    len*sizeof(T), 
                    stack_, 
                    nleft_          ) ) {

      T* result = reinterpret_cast<T*>(stack_);
      stack_ = (char*)stack_ + len*sizeof(T);
      nleft_ -= std::distance( old_stack, 
                               (char*)stack_ );
      return result;

    }

    //throw std::bad_alloc();
    throw std::runtime_error("std::bad_alloc " + msg);

  }

  template <typename T>
  T* aligned_alloc( size_t len, std::string msg ) {
    return aligned_alloc<T>( len, alignof(T), msg );
  }


  inline void* stack() const {return stack_;}
  inline size_t nleft() const { return nleft_; }

};


}
