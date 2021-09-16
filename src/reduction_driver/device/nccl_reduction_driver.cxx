#include "nccl_reduction_driver.hpp"
#include <cstring>
#include <memory>
#include <map>
#include <iostream>

#include <gauxc/util/cuda_util.hpp>
#include "device/type_erased_queue.hpp"

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

cudaStream_t get_cuda_stream_from_optional_args( std::any& args ) {
  cudaStream_t stream = 0;

  if( args.has_value() ) {
    if( auto ptr = std::any_cast<type_erased_queue>( &args ) ) 
    if( auto passed_stream = ptr->queue_as_ptr<util::cuda_stream>() ) {
      stream = *passed_stream;
    }
  }

  return stream;
}


NCCLReductionDriver::NCCLReductionDriver(GAUXC_MPI_CODE(MPI_Comm comm)) :
  DeviceReductionDriver(GAUXC_MPI_CODE(comm)),
  nccl_comm_( std::make_shared<util::nccl_comm>(comm) ){ }


NCCLReductionDriver::~NCCLReductionDriver() noexcept = default;
NCCLReductionDriver::NCCLReductionDriver(const NCCLReductionDriver&) = default;


void NCCLReductionDriver::allreduce_typeerased( const void* src, void* dest, 
  size_t size, ReductionOp op, std::type_index idx, std::any optional_args )  {

  auto stream = get_cuda_stream_from_optional_args( optional_args );

  auto synchronize = [&]() {
    if( stream == 0 ) cudaDeviceSynchronize();
    else              cudaStreamSynchronize(stream);
  };

  synchronize();
  auto err = ncclAllReduce( src, dest, size, get_nccl_datatype(idx), 
    get_nccl_op(op), *nccl_comm_, 0 );

  if( err != ncclSuccess ) throw std::runtime_error("NCCL FAILED");
  synchronize();

}
void NCCLReductionDriver::allreduce_inplace_typeerased( void* data, size_t size,
  ReductionOp op, std::type_index idx, std::any optional_args) {

  auto stream = get_cuda_stream_from_optional_args( optional_args );

  auto synchronize = [&]() {
    if( stream == 0 ) cudaDeviceSynchronize();
    else              cudaStreamSynchronize(stream);
  };

  synchronize();
  auto err = ncclAllReduce( data, data, size, get_nccl_datatype(idx),
    get_nccl_op(op), *nccl_comm_, stream );

  if( err != ncclSuccess ) throw std::runtime_error("NCCL FAILED");
  synchronize();


}

std::unique_ptr<detail::ReductionDriverImpl> NCCLReductionDriver::clone() {
  return std::make_unique<NCCLReductionDriver>(*this);
}


}

