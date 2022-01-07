#include "reduction_driver_impl.hpp"

namespace GauXC::detail {


ReductionDriverImpl::ReductionDriverImpl( GAUXC_MPI_CODE(MPI_Comm comm) ) 
  GAUXC_MPI_CODE( : comm_(comm) ) {}

ReductionDriverImpl::~ReductionDriverImpl() noexcept = default;
ReductionDriverImpl::ReductionDriverImpl(const ReductionDriverImpl& ) = default;

#ifdef GAUXC_ENABLE_MPI
MPI_Comm ReductionDriverImpl::comm() const {
  return comm_;
}
#endif

}
