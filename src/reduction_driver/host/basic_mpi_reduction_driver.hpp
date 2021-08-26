#pragma once
#include "host_reduction_driver.hpp"

namespace GauXC {

struct BasicMPIReductionDriver : public HostReductionDriver {

  BasicMPIReductionDriver(GAUXC_MPI_CODE(MPI_Comm comm));
  virtual ~BasicMPIReductionDriver() noexcept;
  BasicMPIReductionDriver(const BasicMPIReductionDriver& );

  void allreduce_typeerased( const void*, void*, size_t, ReductionOp, std::type_index ) override;
  void allreduce_inplace_typeerased( void*, size_t, ReductionOp, std::type_index ) override;
  
  std::unique_ptr<detail::ReductionDriverImpl> clone() override;

};

}
