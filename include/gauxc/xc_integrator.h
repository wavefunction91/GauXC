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

#ifdef __cplusplus
#include <cstdint>
#include <cstdbool>
#else
#include <stdint.h>
#include <stdbool.h>
#endif
#include <gauxc/types.h>
#include <gauxc/status.h>
#include <gauxc/functional.h>
#include <gauxc/load_balancer.h>

#ifdef __cplusplus
extern "C" {
namespace GauXC::C {
#endif

/**
 * @brief GauXC C API XCIntegrator handle.
 */
typedef struct GauXCIntegrator {
  GauXCHeader hdr; ///< Header for internal use.
  void* ptr; ///< Pointer to the XCIntegrator instance.
} GauXCIntegrator;

/**
 * @brief Delete an XCIntegrator instance.
 * @param status Status object to capture any errors.
 * @param integrator Handle to the XCIntegrator to delete.
 */
extern void gauxc_integrator_delete(
  GauXCStatus* status,
  GauXCIntegrator* integrator
);

/**
 * @brief Create a new XCIntegrator instance.
 * @param status Status object to capture any errors.
 * @param functional Handle to the Functional.
 * @param lb Handle to the LoadBalancer.
 * @param execution_space Execution space to use.
 * @param integrator_input_type Name of the integrator input type.
 * @param integrator_kernel_name Name of the integrator kernel.
 * @param local_work_kernel_name Name of the local work kernel.
 * @param reduction_kernel_name Name of the reduction kernel.
 * @return Handle to the created XCIntegrator.
 */
extern GauXCIntegrator gauxc_integrator_new(
  GauXCStatus* status,
  const GauXCFunctional functional,
  const GauXCLoadBalancer lb,
  enum GauXC_ExecutionSpace execution_space,
  const char* integrator_input_type,
  const char* integrator_kernel_name,
  const char* local_work_kernel_name,
  const char* reduction_kernel_name
);

/**
 * @brief Integrate the density matrix to get the number of electrons.
 * @param status Status object to capture any errors.
 * @param integrator Handle to the XCIntegrator.
 * @param m Number of rows in the density matrix.
 * @param n Number of columns in the density matrix.
 * @param density_matrix Pointer to the density matrix data.
 * @param ldp Leading dimension of the density matrix.
 * @param den Pointer to store the number of electrons.
 */
extern void gauxc_integrator_integrate_den(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  const int64_t m,
  const int64_t n,
  const double* density_matrix,
  const int64_t ldp,
  double* den
);

/**
 * @brief Evaluate the exchange-correlation energy for RKS.
 * @param status Status object to capture any errors.
 * @param integrator Handle to the XCIntegrator.
 * @param m Number of rows in the density matrix.
 * @param n Number of columns in the density matrix.
 * @param density_matrix Pointer to the density matrix data.
 * @param ldp Leading dimension of the density matrix.
 * @param exc Pointer to store the exchange-correlation energy.
 */
extern void gauxc_integrator_eval_exc_rks(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  const int64_t m,
  const int64_t n,
  const double* density_matrix,
  const int64_t ldp,
  double* exc
);

/**
 * @brief Evaluate the exchange-correlation energy for UKS.
 * @param status Status object to capture any errors.
 * @param integrator Handle to the XCIntegrator.
 * @param m Number of rows in the density matrix.
 * @param n Number of columns in the density matrix.
 * @param density_matrix_s Pointer to the density matrix data for total density.
 * @param ldp_s Leading dimension of the total density matrix.
 * @param density_matrix_z Pointer to the density matrix data for spin density.
 * @param ldp_z Leading dimension of the spin density matrix.
 * @param exc Pointer to store the exchange-correlation energy.
 */
extern void gauxc_integrator_eval_exc_uks(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  const int64_t m,
  const int64_t n,
  const double* density_matrix_s,
  const int64_t ldp_s,
  const double* density_matrix_z,
  const int64_t ldp_z,
  double* exc
);

/**
 * @brief Evaluate the exchange-correlation energy for GKS.
 * @param status Status object to capture any errors.
 * @param integrator Handle to the XCIntegrator.
 * @param m Number of rows in the density matrix.
 * @param n Number of columns in the density matrix.
 * @param density_matrix_s Pointer to the density matrix data for total density.
 * @param ldp_s Leading dimension of the total density matrix.
 * @param density_matrix_z Pointer to the density matrix data for spin z density.
 * @param ldp_z Leading dimension of the spin z density matrix.
 * @param density_matrix_y Pointer to the density matrix data for spin y component.
 * @param ldp_y Leading dimension of the spin y density matrix.
 * @param density_matrix_x Pointer to the density matrix data for spin x component.
 * @param ldp_x Leading dimension of the spin x density matrix.
 * @param exc Pointer to store the exchange-correlation energy.
 */
extern void gauxc_integrator_eval_exc_gks(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  const int64_t m,
  const int64_t n,
  const double* density_matrix_s,
  const int64_t ldp_s,
  const double* density_matrix_z,
  const int64_t ldp_z,
  const double* density_matrix_y,
  const int64_t ldp_y,
  const double* density_matrix_x,
  const int64_t ldp_x,
  double* exc
);

/**
 * @brief Evaluate the exchange-correlation energy and potential for RKS.
 * @param status Status object to capture any errors.
 * @param integrator Handle to the XCIntegrator.
 * @param m Number of rows in the density matrix.
 * @param n Number of columns in the density matrix.
 * @param density_matrix Pointer to the density matrix data.
 * @param ldp Leading dimension of the density matrix.
 * @param exc Pointer to store the exchange-correlation energy.
 * @param vxc_matrix Matrix container to store the exchange-correlation potential.
 * @param vxc_ld Leading dimension of the exchange-correlation potential matrix.
 */
extern void gauxc_integrator_eval_exc_vxc_rks(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  const int64_t m,
  const int64_t n,
  const double* density_matrix,
  const int64_t ldp,
  double* exc,
  double* vxc_matrix,
  const int64_t vxc_ld
);

/**
 * @brief Evaluate the exchange-correlation energy and potential for UKS.
 * @param status Status object to capture any errors.
 * @param integrator Handle to the XCIntegrator.
 * @param m Number of rows in the density matrix.
 * @param n Number of columns in the density matrix.
 * @param density_matrix_s Pointer to the density matrix data for total density.
 * @param ldp_s Leading dimension of the total density matrix.
 * @param density_matrix_z Pointer to the density matrix data for spin density.
 * @param ldp_z Leading dimension of the spin density matrix.
 * @param exc Pointer to store the exchange-correlation energy.
 * @param vxc_matrix_s Matrix container to store the exchange-correlation potential for total density.
 * @param vxc_ld_s Leading dimension of the exchange-correlation potential matrix for total density.
 * @param vxc_matrix_z Matrix container to store the exchange-correlation potential for spin density.
 * @param vxc_ld_z Leading dimension of the exchange-correlation potential matrix for spin density.
 */
extern void gauxc_integrator_eval_exc_vxc_uks(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  const int64_t m,
  const int64_t n,
  const double* density_matrix_s,
  const int64_t ldp_s,
  const double* density_matrix_z,
  const int64_t ldp_z,
  double* exc,
  double* vxc_matrix_s,
  const int64_t vxc_ld_s,
  double* vxc_matrix_z,
  const int64_t vxc_ld_z
);

/**
 * @brief Evaluate the exchange-correlation energy and potential for GKS.
 * @param status Status object to capture any errors.
 * @param integrator Handle to the XCIntegrator.
 * @param m Number of rows in the density matrix.
 * @param n Number of columns in the density matrix.
 * @param density_matrix_s Pointer to the density matrix data for total density.
 * @param ldp_s Leading dimension of the total density matrix.
 * @param density_matrix_z Pointer to the density matrix data for spin z density.
 * @param ldp_z Leading dimension of the spin z density matrix.
 * @param density_matrix_y Pointer to the density matrix data for spin y component.
 * @param ldp_y Leading dimension of the spin y density matrix.
 * @param density_matrix_x Pointer to the density matrix data for spin x component.
 * @param ldp_x Leading dimension of the spin x density matrix.
 * @param exc Pointer to store the exchange-correlation energy.
 * @param vxc_matrix_s Matrix container to store the exchange-correlation potential for total density.
 * @param vxc_ld_s Leading dimension of the exchange-correlation potential matrix for total density.
 * @param vxc_matrix_z Matrix container to store the exchange-correlation potential for spin z density.
 * @param vxc_ld_z Leading dimension of the exchange-correlation potential matrix for spin z density.
 * @param vxc_matrix_y Matrix container to store the exchange-correlation potential for spin y component.
 * @param vxc_ld_y Leading dimension of the exchange-correlation potential matrix for spin y component.
 * @param vxc_matrix_x Matrix container to store the exchange-correlation potential for spin x component.
 * @param vxc_ld_x Leading dimension of the exchange-correlation potential matrix for spin x component.
 */
extern void gauxc_integrator_eval_exc_vxc_gks(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  const int64_t m,
  const int64_t n,
  const double* density_matrix_s,
  const int64_t ldp_s,
  const double* density_matrix_z,
  const int64_t ldp_z,
  const double* density_matrix_y,
  const int64_t ldp_y,
  const double* density_matrix_x,
  const int64_t ldp_x,
  double* exc,
  double* vxc_matrix_s,
  const int64_t vxc_ld_s,
  double* vxc_matrix_z,
  const int64_t vxc_ld_z,
  double* vxc_matrix_y,
  const int64_t vxc_ld_y,
  double* vxc_matrix_x,
  const int64_t vxc_ld_x
);

#ifdef __cplusplus
} // namespace GauXC::C
} // extern "C"
#endif