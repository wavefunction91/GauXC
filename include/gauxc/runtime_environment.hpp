#pragma once
#include <gauxc/gauxc_config.hpp>
#include <memory>
#ifdef GAUXC_ENABLE_MPI
#include <mpi.h>
#endif

namespace GauXC {

namespace detail {
  class RuntimeEnvironmentImpl;
}

class RuntimeEnvironment {

protected:

  using pimpl_type = detail::RuntimeEnvironmentImpl;
  using pimpl_ptr_type = std::unique_ptr<pimpl_type>;
  pimpl_ptr_type pimpl_;
  RuntimeEnvironment( pimpl_ptr_type&& ptr );

public:

  explicit RuntimeEnvironment(GAUXC_MPI_CODE(MPI_Comm comm));
  virtual ~RuntimeEnvironment() noexcept;

  RuntimeEnvironment( const RuntimeEnvironment& ) = delete;
  RuntimeEnvironment( RuntimeEnvironment&& ) noexcept;

  GAUXC_MPI_CODE(MPI_Comm comm() const;)
  int comm_rank() const;
  int comm_size() const;

};

#ifdef GAUXC_ENABLE_DEVICE
class DeviceRuntimeEnvironment : public RuntimeEnvironment {

private:
  using parent_type = RuntimeEnvironment;

  using parent_type::pimpl_type;
  using parent_type::pimpl_ptr_type;

public:

  DeviceRuntimeEnvironment(GAUXC_MPI_CODE(MPI_Comm comm,) void* mem, 
    size_t mem_sz);
  DeviceRuntimeEnvironment(GAUXC_MPI_CODE(MPI_Comm));

  ~DeviceRuntimeEnvironment() noexcept;
  DeviceRuntimeEnvironment( const DeviceRuntimeEnvironment& ) = delete;
  DeviceRuntimeEnvironment( DeviceRuntimeEnvironment&& ) noexcept;

  void* device_memory();
  size_t device_memory_size();

};
#endif

}
