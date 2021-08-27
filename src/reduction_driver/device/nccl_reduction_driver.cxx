#include "nccl_reduction_driver.hpp"
#include <cstring>
#include <memory>
#include <map>
#include <iostream>

namespace GauXC {

ncclDataType_t get_nccl_datatype( std::type_index idx ) {

  static std::map<std::type_index, ncclDataType_t> map {
    {std::type_index(typeid(double)), ncclDouble},
    {std::type_index(typeid(float)),  ncclFloat}
  };

  return map.at(idx);

}

ncclRedOp_t get_nccl_op( ReductionOp op ) {

  static std::map< ReductionOp, ncclRedOp_t > map {
    { ReductionOp::Sum, ncclSum }
  };

  return map.at(op);

}


NCCLReductionDriver::NCCLReductionDriver(GAUXC_MPI_CODE(MPI_Comm comm)) :
  DeviceReductionDriver(GAUXC_MPI_CODE(comm)),
  nccl_comm_( std::make_shared<util::nccl_comm>(comm) ){ }


NCCLReductionDriver::~NCCLReductionDriver() noexcept = default;
NCCLReductionDriver::NCCLReductionDriver(const NCCLReductionDriver&) = default;


void NCCLReductionDriver::allreduce_typeerased( const void* src, void* dest, 
  size_t size, ReductionOp op, std::type_index idx, std::any )  {

  cudaDeviceSynchronize();
  auto err = ncclAllReduce( src, dest, size, get_nccl_datatype(idx), 
    get_nccl_op(op), *nccl_comm_, 0 );
  cudaDeviceSynchronize();

  if( err != ncclSuccess ) throw std::runtime_error("NCCL FAILED");

}
void NCCLReductionDriver::allreduce_inplace_typeerased( void* data, size_t size,
  ReductionOp op, std::type_index idx, std::any) {

  cudaDeviceSynchronize();
  auto err =ncclAllReduce( data, data, size, get_nccl_datatype(idx),
    get_nccl_op(op), *nccl_comm_, 0 );
  cudaDeviceSynchronize();

  if( err != ncclSuccess ) throw std::runtime_error("NCCL FAILED");

}

std::unique_ptr<detail::ReductionDriverImpl> NCCLReductionDriver::clone() {
  return std::make_unique<NCCLReductionDriver>(*this);
}


}

