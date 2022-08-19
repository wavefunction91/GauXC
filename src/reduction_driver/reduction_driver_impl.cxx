#include "reduction_driver_impl.hpp"

namespace GauXC::detail {


ReductionDriverImpl::ReductionDriverImpl( const RuntimeEnvironment& rt ) 
  : runtime_(rt){}

ReductionDriverImpl::~ReductionDriverImpl() noexcept = default;
ReductionDriverImpl::ReductionDriverImpl(const ReductionDriverImpl& ) = default;

}
