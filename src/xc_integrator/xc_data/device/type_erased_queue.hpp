#pragma once
#include <any>
#include <memory>
#include <type_traits>

namespace GauXC {
//namespace detail {
//  template <typename T>
//  struct is_shared_ptr : public std::false_type {};
//
//  template <typename T>
//  struct is_shared_ptr<std::shared_ptr<T>> : public std::true_type {};
//
//
//  template <typename T, typename U = void>
//  struct enable_if_not_shared_ptr : public std::enable_if< not is_shared_ptr<T>::value, U> {};
//
//}

class type_erased_queue {

  std::any queue_;

public:

  type_erased_queue() = default;

  template <typename T>
  type_erased_queue( std::shared_ptr<T> q ) : queue_( std::move(q) ) { }

  //template <typename T, typename = typename detail::enable_if_not_shared_ptr<std::decay_t<T>>::type >
  //type_erased_queue( T&& q ) :
  //  type_erased_queue(std::make_shared<std::decay_t<T>>(std::forward<T>(q))) { }

  template <typename T>
  inline const T* queue_as_ptr() const { 
    if( !queue_.has_value() ) return nullptr;
    if( auto q_ptr = std::any_cast< std::shared_ptr<T> >( &queue_ ) ) {
      return q_ptr->get();
    } else {
      return nullptr;
    }
  }

  template <typename T>
  inline T* queue_as_ptr() { 
    if( auto q_ptr = std::any_cast< std::shared_ptr<T> >( &queue_ ) ) {
      return q_ptr->get();
    } else {
      return nullptr;
    }
  }

  template <typename T>
  inline const T& queue_as() const {
    auto ptr = queue_as_ptr<T>();
    if( not ptr ) throw std::bad_any_cast();
    return *ptr;
  }

  template <typename T>
  inline T& queue_as() {
    auto ptr = queue_as_ptr<T>();
    if( not ptr ) throw std::bad_any_cast();
    return *ptr;
  }

};

}
