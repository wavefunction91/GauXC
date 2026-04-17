/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy).
 *
 * (c) 2024-2025, Microsoft Corporation
 *
 * All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include <gauxc/runtime_environment/fwd.hpp>
#include <memory>
#include <gauxc/util/mpi.hpp>

namespace GauXC {

namespace detail {
  /// Forward declaration of RuntimeEnvironment implementation.
  class RuntimeEnvironmentImpl;
  #ifdef GAUXC_HAS_DEVICE
  /// Convert a RuntimeEnvironment to a DeviceRuntimeEnvironment.
  DeviceRuntimeEnvironment as_device_runtime( const RuntimeEnvironment& );
  #endif
}

/**
 *  @brief Runtime environment for GauXC computations.
 *
 *  Encapsulates communication context (MPI) and provides process-local
 *  information such as rank and size. Serves as the base class for
 *  device-specific runtime environments.
 */
class RuntimeEnvironment {

protected:

#ifdef GAUXC_HAS_DEVICE
  friend DeviceRuntimeEnvironment 
    detail::as_device_runtime(const RuntimeEnvironment&); 
#endif

  using pimpl_type = detail::RuntimeEnvironmentImpl;       ///< Implementation type.
  using pimpl_ptr_type = std::shared_ptr<pimpl_type>;      ///< Shared pointer to impl.
  pimpl_ptr_type pimpl_;                                   ///< Handle for implementation object.

  /// Construct from existing implementation pointer.
  RuntimeEnvironment( pimpl_ptr_type ptr );

public:

  /**
   *  @brief Construct a runtime environment.
   *  @param comm MPI communicator (if MPI is enabled).
   */
  explicit RuntimeEnvironment(GAUXC_MPI_CODE(MPI_Comm comm));

  /// Destructor.
  virtual ~RuntimeEnvironment() noexcept;

  /// Copy constructor.
  RuntimeEnvironment( const RuntimeEnvironment& );
  /// Move constructor.
  RuntimeEnvironment( RuntimeEnvironment&& ) noexcept;

  /// Get the MPI communicator (if MPI is enabled).
  GAUXC_MPI_CODE(MPI_Comm comm() const;)

  /// Get the rank of this process in the communicator.
  int comm_rank() const;

  /// Get the total number of processes in the communicator.
  int comm_size() const;

  /// Get the shared usage count of the implementation.
  int shared_usage_count() const;

};

#ifdef GAUXC_HAS_DEVICE
/**
 *  @brief Device-enabled runtime environment for GPU computations.
 *
 *  Extends RuntimeEnvironment with device memory management and backend access
 *  for GPU-accelerated XC integration.
 */
class DeviceRuntimeEnvironment : public RuntimeEnvironment {

private:

  using parent_type = RuntimeEnvironment;
  friend DeviceRuntimeEnvironment 
    detail::as_device_runtime(const RuntimeEnvironment&); 

  using parent_type::pimpl_type;
  using parent_type::pimpl_ptr_type;

  /// Construct from existing implementation pointer.
  DeviceRuntimeEnvironment( pimpl_ptr_type ptr );

public:

  /**
   *  @brief Construct with explicit device memory buffer.
   *  @param comm MPI communicator (if MPI is enabled).
   *  @param mem  Pointer to device memory buffer.
   *  @param mem_sz Size of device memory buffer in bytes.
   */
  DeviceRuntimeEnvironment(GAUXC_MPI_CODE(MPI_Comm comm,) void* mem, 
    size_t mem_sz);

  /**
   *  @brief Construct with automatic device memory allocation.
   *  @param comm MPI communicator (if MPI is enabled).
   *  @param fill_fraction Fraction of available device memory to allocate.
   */
  DeviceRuntimeEnvironment(GAUXC_MPI_CODE(MPI_Comm,) double fill_fraction);

  /// Destructor.
  ~DeviceRuntimeEnvironment() noexcept;

  /// Copy constructor.
  DeviceRuntimeEnvironment( const DeviceRuntimeEnvironment& );
  /// Move constructor.
  DeviceRuntimeEnvironment( DeviceRuntimeEnvironment&& ) noexcept;

  /// Get pointer to device memory buffer.
  void* device_memory() const ;

  /// Get size of device memory buffer in bytes.
  size_t device_memory_size() const ;

  /// Check if this environment owns the device memory.
  bool owns_memory() const;

  /// Get the device backend instance.
  DeviceBackend* device_backend() const;

  /// Release the device memory buffer.
  void release_buffer();

  /**
   *  @brief Set an external device memory buffer.
   *  @param m  Pointer to device memory.
   *  @param sz Size of device memory buffer in bytes.
   */
  void set_buffer(void* m, size_t sz);
};
#endif

}
