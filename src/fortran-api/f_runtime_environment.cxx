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
#include <gauxc/runtime_environment.h>
#include <gauxc/util/mpi.h>

namespace GauXC::C {
extern "C" {

GauXCRuntimeEnvironment gauxc_runtime_environment_new_f( 
  GauXCStatus* status
  GAUXC_MPI_CODE(, MPI_Comm comm) 
) {
  return gauxc_runtime_environment_new(status GAUXC_MPI_CODE(, MPI_Comm_f2c(comm)));
}

#ifdef GAUXC_HAS_DEVICE
GauXCRuntimeEnvironment gauxc_device_runtime_environment_new_f( 
  GauXCStatus* status,
  GAUXC_MPI_CODE(MPI_Comm comm,)
  double fill_fraction
) {
  return gauxc_device_runtime_environment_new(status GAUXC_MPI_CODE(, MPI_Comm_f2c(comm)), fill_fraction);
}

GauXCRuntimeEnvironment gauxc_device_runtime_environment_new_mem_f( 
  GauXCStatus* status,
  GAUXC_MPI_CODE(MPI_Comm comm,)
  void* mem,
  size_t mem_sz
) {
  return gauxc_device_runtime_environment_new_mem(status GAUXC_MPI_CODE(, MPI_Comm_f2c(comm)), mem, mem_sz);
}
#endif

} // extern "C"
} // namespace GauXC::C
