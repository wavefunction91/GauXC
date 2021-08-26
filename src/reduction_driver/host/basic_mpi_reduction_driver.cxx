#include "basic_mpi_reduction_driver.hpp"
#include <cstring>
#include <memory>
#include <map>

namespace GauXC {

MPI_Datatype get_datatype( std::type_index idx ) {

  static std::map<std::type_index, MPI_Datatype> map {
    {std::type_index(typeid(double)), MPI_DOUBLE},
    {std::type_index(typeid(float)),  MPI_FLOAT}
  };

  return map.at(idx);

}

MPI_Op get_op( ReductionOp op ) {

  static std::map< ReductionOp, MPI_Op > map {
    { ReductionOp::Sum, MPI_SUM }
  };

  return map.at(op);

}


BasicMPIReductionDriver::BasicMPIReductionDriver(GAUXC_MPI_CODE(MPI_Comm comm)) :
  HostReductionDriver(GAUXC_MPI_CODE(comm)) { }


BasicMPIReductionDriver::~BasicMPIReductionDriver() noexcept = default;
BasicMPIReductionDriver::BasicMPIReductionDriver(const BasicMPIReductionDriver&) = default;


void BasicMPIReductionDriver::allreduce_typeerased( const void* src, void* dest, 
  size_t size, ReductionOp op, std::type_index idx )  {

#ifdef GAUXC_ENABLE_MPI
  int world_size;
  MPI_Comm_size( comm_, &world_size );
#else
  int world_size = 1;
#endif

  if( world_size == 1 ) {
    std::memcpy( dest, src, size ); 
  } else  {
    #ifdef GAUXC_ENABLE_MPI 
    MPI_Allreduce( src, dest, size, get_datatype(idx), get_op(op), comm_ );
    #endif
  }


}
void BasicMPIReductionDriver::allreduce_inplace_typeerased( void* data, size_t size,
  ReductionOp op, std::type_index idx) {

#ifdef GAUXC_ENABLE_MPI
  int world_size;
  MPI_Comm_size( comm_, &world_size );
#else
  int world_size = 1;
#endif

  if(world_size > 1) {
    #ifdef GAUXC_ENABLE_MPI
    // Test of communicator is an inter-communicator
    int inter_flag;
    MPI_Comm_test_inter( comm_, &inter_flag );

    // Reduce in place
    if( not inter_flag ) {
      MPI_Allreduce( MPI_IN_PLACE, data, size, get_datatype(idx), get_op(op), comm_ );

    // Cannot reduce in place
    } else {
      std::allocator<std::byte> alloc;
      auto* tmp = alloc.allocate( size );
      std::memcpy(tmp, data, size);
      allreduce_typeerased( tmp, data, size, op, idx );
      alloc.deallocate( tmp, size );
    }
    #endif
  }
}

std::unique_ptr<detail::ReductionDriverImpl> BasicMPIReductionDriver::clone() {
  return std::make_unique<BasicMPIReductionDriver>(*this);
}


}

