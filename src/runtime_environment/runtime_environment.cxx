#include <gauxc/runtime_environment.hpp>
#include "runtime_environment_impl.hpp"

namespace GauXC {

RuntimeEnvironment::RuntimeEnvironment( pimpl_ptr_type&& ptr ) :
  pimpl_(std::move(ptr)) {}

RuntimeEnvironment::RuntimeEnvironment(GAUXC_MPI_CODE(MPI_Comm c)) :
  RuntimeEnvironment( std::make_unique<detail::RuntimeEnvironmentImpl>(GAUXC_MPI_CODE(c)) ) {}

RuntimeEnvironment::~RuntimeEnvironment() noexcept = default;

#ifdef GAUXC_ENABLE_MPI
MPI_Comm RuntimeEnvironment::comm() const {
  return pimpl_->comm();
}
#endif

int RuntimeEnvironment::comm_rank() const {
  return pimpl_->comm_rank();
}

int RuntimeEnvironment::comm_size() const {
  return pimpl_->comm_size();
}

}
