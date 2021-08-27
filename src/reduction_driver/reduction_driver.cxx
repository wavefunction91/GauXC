#include "reduction_driver_impl.hpp"

namespace GauXC {

ReductionDriver::ReductionDriver( std::unique_ptr<pimpl_type>&& pimpl ): 
  pimpl_( std::move(pimpl) ) { }

ReductionDriver::ReductionDriver() : ReductionDriver( nullptr ) { }

ReductionDriver::ReductionDriver( const ReductionDriver& other ) :
  ReductionDriver(other.pimpl_->clone()){ }

ReductionDriver::ReductionDriver( ReductionDriver&& ) noexcept = default;
              
ReductionDriver::~ReductionDriver() noexcept = default;


void ReductionDriver::allreduce_typeerased( const void* src, void* dest, size_t size, ReductionOp op, std::type_index idx, std::any optional_args ) {
  if( not pimpl_ ) throw std::runtime_error("Not Initialized");
  pimpl_->allreduce_typeerased(src, dest, size, op, idx, optional_args);
}

void ReductionDriver::allreduce_inplace_typeerased( void* data, size_t size, ReductionOp op, std::type_index idx, std::any optional_args ) {
  if( not pimpl_ ) throw std::runtime_error("Not Initialized");
  pimpl_->allreduce_inplace_typeerased(data, size, op, idx, optional_args);
}

bool ReductionDriver::takes_host_memory() const {
  if( not pimpl_ ) throw std::runtime_error("Not Initialized");
  return pimpl_->takes_host_memory();
}
bool ReductionDriver::takes_device_memory() const {
  if( not pimpl_ ) throw std::runtime_error("Not Initialized");
  return pimpl_->takes_device_memory();
}


#ifdef GAUXC_ENABLE_MPI
MPI_Comm ReductionDriver::comm() const {
  if( not pimpl_ ) throw std::runtime_error("Not Initialized");
  return pimpl_->comm();
}
#endif


}
