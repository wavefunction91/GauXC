#include "device_runtime_environment_impl.hpp"
#include <gauxc/exceptions.hpp>

namespace GauXC {

DeviceRuntimeEnvironment::DeviceRuntimeEnvironment(
  GAUXC_MPI_CODE(MPI_Comm c,) void* p, size_t sz ) :
  RuntimeEnvironment(
    detail::make_device_runtime( GAUXC_MPI_CODE(c,) p,sz)
  ) {}

DeviceRuntimeEnvironment::DeviceRuntimeEnvironment(
  GAUXC_MPI_CODE(MPI_Comm c), double ff) :
  RuntimeEnvironment(detail::make_device_runtime(GAUXC_MPI_CODE(c,)ff)) {}

DeviceRuntimeEnvironment::~DeviceRuntimeEnvironment() noexcept = default;

auto* device_runtime_pimpl_cast(detail::RuntimeEnvironmentImpl* ptr) {
  auto dp = dynamic_cast<detail::DeviceRuntimeEnvironmentImpl*>(ptr);
  if(!dp) GAUXC_GENERIC_EXCEPTION("Not A Device Implemention");
  return dp;
}

void* DeviceRuntimeEnvironment::device_memory() const {
  return device_runtime_pimpl_cast(pimpl_.get())->device_memory();
}
size_t DeviceRuntimeEnvironment::device_memory_size() const {
  return device_runtime_pimpl_cast(pimpl_.get())->device_memory_size();
}
DeviceBackend* DeviceRuntimeEnvironment::device_backend() const {
  return device_runtime_pimpl_cast(pimpl_.get())->device_backend();
}


}
