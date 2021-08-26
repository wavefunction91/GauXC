#pragma once
#include <gauxc/reduction_driver.hpp>

namespace GauXC  {
namespace detail {

class ReductionDriverImpl {

protected: 

  GAUXC_MPI_CODE( MPI_Comm comm_; )

public:

  #ifdef GAUXC_ENABLE_MPI
  ReductionDriverImpl() = delete;
  #endif

  ReductionDriverImpl( GAUXC_MPI_CODE(MPI_Comm comm) );
  virtual ~ReductionDriverImpl() noexcept;
  ReductionDriverImpl( const ReductionDriverImpl& );

  virtual void allreduce_typeerased( const void*, void*, size_t, ReductionOp, std::type_index ) = 0;
  virtual void allreduce_inplace_typeerased( void*, size_t, ReductionOp, std::type_index ) = 0;

  virtual bool takes_host_memory() const = 0;
  virtual bool takes_device_memory() const = 0;

#ifdef GAUXC_ENABLE_MPI
  MPI_Comm comm() const;
#endif

  virtual std::unique_ptr<ReductionDriverImpl> clone() = 0;
};

}
}
