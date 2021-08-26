#include "reduction_driver_impl.hpp"
#include "host/basic_mpi_reduction_driver.hpp"
#include <algorithm>
#include <iostream>

namespace GauXC {

std::shared_ptr<ReductionDriver> ReductionDriverFactory::get_shared_instance(
  GAUXC_MPI_CODE(MPI_Comm comm,) std::string kernel_name ) {

  std::transform(kernel_name.begin(), kernel_name.end(), 
    kernel_name.begin(), ::toupper );

  std::unique_ptr<detail::ReductionDriverImpl> ptr = nullptr;

  if( kernel_name == "DEFAULT" ) kernel_name = "BASICMPI";

  if( kernel_name == "BASICMPI" )
    ptr = std::make_unique<BasicMPIReductionDriver>(GAUXC_MPI_CODE(comm));

  if( !ptr ) throw std::runtime_error("Unknown Reduction Driver");

  return std::make_shared<ReductionDriver>(std::move(ptr));

}

}
