#pragma once
#include <gauxc/gauxc_config.hpp>
#include <memory>

namespace GauXC  {
namespace detail {

template <typename T>
std::shared_ptr<std::decay_t<T>> forward_as_shared_ptr( const T& t ) {
  return std::make_shared<std::decay_t<T>>( t );
}

//template <typename T>
//std::shared_ptr<std::decay_t<T>> forward_as_shared_ptr( T& t ) {
//  std::cout << "Resolving Ref Copy Forward" << std::endl;
//  return std::make_shared<std::decay_t<T>>( t );
//}
//
//template <typename T>
//std::shared_ptr<std::decay_t<T>> forward_as_shared_ptr( T&& t ) {
//  std::cout << "Resolving Move Forward" << std::endl;
//  return std::make_shared<std::decay_t<T>>( std::move(t) );
//}

template <typename T>
std::shared_ptr<T> forward_as_shared_ptr( std::shared_ptr<T> ptr ) {
  return ptr;
}

// Disable forward for MPI_Comm
#ifdef GAUXC_ENABLE_MPI
MPI_Comm forward_as_shared_ptr( MPI_Comm comm ) {
  return comm;
}
#endif

}
}
