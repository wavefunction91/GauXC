#pragma once
#include <gauxc/reduction_driver.hpp>

namespace GauXC  {
namespace detail {

class ReductionDriverImpl {

protected: 

  const RuntimeEnvironment& runtime_;

public:

  ReductionDriverImpl() = delete;
  ReductionDriverImpl( const RuntimeEnvironment& rt);

  virtual ~ReductionDriverImpl() noexcept;
  ReductionDriverImpl( const ReductionDriverImpl& );

  virtual void allreduce_typeerased( const void*, void*, size_t, ReductionOp, std::type_index, std::any ) = 0;
  virtual void allreduce_inplace_typeerased( void*, size_t, ReductionOp, std::type_index, std::any ) = 0;

  virtual bool takes_host_memory() const = 0;
  virtual bool takes_device_memory() const = 0;

  virtual std::unique_ptr<ReductionDriverImpl> clone() = 0;
};

}
}
