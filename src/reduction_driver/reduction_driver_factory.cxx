#include "reduction_driver_impl.hpp"
#include "host/basic_mpi_reduction_driver.hpp"

#ifdef GAUXC_ENABLE_NCCL
#include "device/nccl_reduction_driver.hpp"
#endif


#include <algorithm>
#include <iostream>
#include <gauxc/exceptions.hpp>

namespace GauXC {

std::shared_ptr<ReductionDriver> ReductionDriverFactory::get_shared_instance(
  GAUXC_MPI_CODE(MPI_Comm comm,) std::string kernel_name ) {

  std::transform(kernel_name.begin(), kernel_name.end(), 
    kernel_name.begin(), ::toupper );

  std::unique_ptr<detail::ReductionDriverImpl> ptr = nullptr;

  if( kernel_name == "DEFAULT" ) kernel_name = "BASICMPI";

  if( kernel_name == "BASICMPI" )
    ptr = std::make_unique<BasicMPIReductionDriver>(GAUXC_MPI_CODE(comm));

  #ifdef GAUXC_ENABLE_NCCL
    if( kernel_name == "NCCL" )
      ptr = std::make_unique<NCCLReductionDriver>(GAUXC_MPI_CODE(comm));
  #endif

  if( !ptr ) GAUXC_GENERIC_EXCEPTION("Unknown Reduction Driver " + kernel_name);

  return std::make_shared<ReductionDriver>(std::move(ptr));

}

}
