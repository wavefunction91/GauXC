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
#include "device_reduction_driver.hpp"
#include <oneapi/ccl.hpp>
#include <memory>

namespace GauXC {

namespace util {
struct oneccl_comm {

  std::shared_ptr<ccl::kvs> kvs;
  ccl::communicator         comm;

  inline oneccl_comm( MPI_Comm mpi_comm, ::sycl::device& device,
    ::sycl::context& context ) :
    comm( make_comm(mpi_comm, device, context, kvs) ) { }

  // oneCCL communicators are move-only (they wrap a PIMPL, cf. ccl::communicator
  // in communicator.hpp), matching the NCCL driver's move-only nccl_comm
  oneccl_comm( const oneccl_comm& ) = delete;
  oneccl_comm( oneccl_comm&& ) noexcept = default;

  inline operator ccl::communicator&() { return comm; }

private:

  // Builds the communicator in a helper so the member initializer list above
  // can construct `comm` directly (ccl::communicator has no default ctor).
  // ccl::create_communicator's DeviceType template parameter is deduced as a
  // non-const reference, so device/context are taken by non-const reference
  // here to match the SYCLBackend members (device, context) passed in below.
  static inline ccl::communicator make_comm( MPI_Comm mpi_comm,
    ::sycl::device& device, ::sycl::context& context,
    std::shared_ptr<ccl::kvs>& kvs_out ) {

    int32_t world_rank, world_size;
    MPI_Comm_rank( mpi_comm, &world_rank );
    MPI_Comm_size( mpi_comm, &world_size );

    ccl::kvs::address_type addr;
    if( world_rank == 0 ) {
      kvs_out = ccl::create_main_kvs();
      addr = kvs_out->get_address();
    }
    MPI_Bcast( addr.data(), addr.size(), MPI_BYTE, 0, mpi_comm );
    if( world_rank != 0 ) kvs_out = ccl::create_kvs( addr );

    return ccl::create_communicator( world_size, world_rank, device, context,
      kvs_out );
  }

};
}


struct OneCCLReductionDriver : public DeviceReductionDriver {

  std::shared_ptr<util::oneccl_comm> oneccl_comm_;

  OneCCLReductionDriver(const RuntimeEnvironment& rt);
  virtual ~OneCCLReductionDriver() noexcept;
  OneCCLReductionDriver(const OneCCLReductionDriver& );

  void allreduce_typeerased( const void*, void*, size_t, ReductionOp, std::type_index, std::any) override;
  void allreduce_inplace_typeerased( void*, size_t, ReductionOp, std::type_index, std::any ) override;

  std::unique_ptr<detail::ReductionDriverImpl> clone() override;

};

}
