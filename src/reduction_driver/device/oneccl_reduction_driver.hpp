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
#include <unordered_map>

namespace GauXC {

namespace util {
struct oneccl_comm {

  std::shared_ptr<ccl::kvs> kvs;
  ccl::communicator         comm;

  inline oneccl_comm( MPI_Comm mpi_comm, ::sycl::device& device,
    ::sycl::context& context ) :
    comm( make_comm(mpi_comm, device, context, kvs) ) { }

  // oneCCL communicators are move-only (they wrap a PIMPL, cf. ccl::communicator
  // in communicator.hpp)
  oneccl_comm( const oneccl_comm& ) = delete;
  oneccl_comm( oneccl_comm&& ) noexcept = default;

  inline operator ccl::communicator&() { return comm; }

private:

  // Builds the communicator in a helper so the member initializer list above
  // can construct `comm` directly (ccl::communicator has no default ctor).
  static inline ccl::communicator make_comm( MPI_Comm mpi_comm,
    const ::sycl::device& device, const ::sycl::context& context,
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

    // create_communicator wraps the device in ccl::pair_class<int, ccl::device>
    // internally, so it needs oneCCL's own device/context types. DeviceType is
    // deduced from a non-const lvalue reference, hence the named local.
    ccl::device ccl_dev = ccl::create_device( ::sycl::device(device) );
    return ccl::create_communicator( world_size, world_rank, ccl_dev,
      ccl::create_context( ::sycl::context(context) ), kvs_out );
  }

};
}


struct OneCCLReductionDriver : public DeviceReductionDriver {

  std::shared_ptr<util::oneccl_comm> oneccl_comm_;

  // ccl::stream wrappers around the SYCL queues this driver has been handed,
  // keyed by queue. Cached so that a collective does not rebuild one per call;
  // shared so that copies of the driver (clone()) reuse the same wrappers.
  std::shared_ptr<std::unordered_map<::sycl::queue, ccl::stream>> ccl_streams_;

  OneCCLReductionDriver(const RuntimeEnvironment& rt);
  virtual ~OneCCLReductionDriver() noexcept;
  OneCCLReductionDriver(const OneCCLReductionDriver& );

  ccl::stream& ccl_stream_for( ::sycl::queue& q );

  void allreduce_typeerased( const void*, void*, size_t, ReductionOp, std::type_index, std::any) override;
  void allreduce_inplace_typeerased( void*, size_t, ReductionOp, std::type_index, std::any ) override;

  std::unique_ptr<detail::ReductionDriverImpl> clone() override;

};

}
