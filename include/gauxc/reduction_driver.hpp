#pragma once
#include <memory>
#include <gauxc/gauxc_config.hpp>
#include <typeindex>
#include <any>
#ifdef GAUXC_ENABLE_MPI
#include <mpi.h>
#endif

namespace GauXC {

namespace detail {
  class ReductionDriverImpl;
}

enum class ReductionOp : int {
  Sum
};

class ReductionDriver {

  using pimpl_type = detail::ReductionDriverImpl;
  std::unique_ptr<pimpl_type> pimpl_;

public:

  ReductionDriver();

  ReductionDriver( std::unique_ptr<pimpl_type>&& pimpl );

  ReductionDriver( const ReductionDriver& );
  ReductionDriver( ReductionDriver&& ) noexcept;

  ~ReductionDriver() noexcept;

  template <typename T>
  inline void allreduce( const T* src, T* dest, size_t size, ReductionOp op, std::any optional_args = std::any()) {
    allreduce_typeerased( src, dest, size, op, std::type_index(typeid(T)), optional_args );
  }

  template <typename T>
  inline void allreduce_inplace( T* data, size_t size, ReductionOp op, std::any optional_args = std::any() ) {
    allreduce_inplace_typeerased( data, size, op, std::type_index(typeid(T)), optional_args );
  }

  void allreduce_typeerased( const void*, void*, size_t, ReductionOp, std::type_index, std::any );
  void allreduce_inplace_typeerased( void*, size_t, ReductionOp, std::type_index, std::any );

  bool takes_host_memory() const;
  bool takes_device_memory() const;

#ifdef GAUXC_ENABLE_MPI
  MPI_Comm comm() const;
#endif

};


struct ReductionDriverFactory {
  static ReductionDriver get_instance( 
    GAUXC_MPI_CODE(MPI_Comm comm,) std::string kernel_name );
  static std::shared_ptr<ReductionDriver> get_shared_instance( 
    GAUXC_MPI_CODE(MPI_Comm comm,) std::string kernel_name );
};

}
