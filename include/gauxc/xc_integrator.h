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
#include <gauxc/matrix.h>

#ifdef __cplusplus
namespace GauXC::C {
extern "C" {
#endif

/**
 * @brief GauXC C API XCIntegrator handle.
 */
typedef struct GauXCIntegrator {
  GauXCHeader hdr; ///< Header for internal use.
  void* ptr; ///< Pointer to the XCIntegrator instance.
  bool owned; ///< Whether this instance owns the XCIntegrator.
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
 * @brief GauXC C API XCIntegratorFactory handle.
 */
typedef struct GauXCIntegratorFactory {
  GauXCHeader hdr; ///< Header for internal use.
  void* ptr; ///< Pointer to the XCIntegratorFactory instance.
} GauXCIntegratorFactory;

/**
 * @brief Create a new XCIntegratorFactory instance.
 * @param status Status object to capture any errors.
 * @param execution_space Execution space to use.
 * @param integrator_input_type Name of the integrator input type.
 * @param integrator_kernel_name Name of the integrator kernel.
 * @param local_work_kernel_name Name of the local work kernel.
 * @param reduction_kernel_name Name of the reduction kernel.
 * @return Handle to the created XCIntegratorFactory.
 */
extern GauXCIntegratorFactory gauxc_integrator_factory_new(
    GauXCStatus* status,
    enum GauXC_ExecutionSpace execution_space,
    const char* integrator_input_type,
    const char* integrator_kernel_name,
    const char* local_work_kernel_name,
    const char* reduction_kernel_name
);

/**
 * @brief Delete an XCIntegratorFactory instance.
 * @param status Status object to capture any errors.
 * @param integrator_factory Handle to the XCIntegratorFactory to delete.
 */
extern void gauxc_integrator_factory_delete(
    GauXCStatus* status,
    GauXCIntegratorFactory* integrator_factory
);

/**
 * @brief Get an XCIntegrator instance from an XCIntegratorFactory.
 * @param status Status object to capture any errors.
 * @param integrator_factory Handle to the XCIntegratorFactory.
 * @param func Handle to the XCFunctional.
 * @param lb Handle to the LoadBalancer.
 * @return Handle to the created XCIntegrator.
 */
extern GauXCIntegrator gauxc_integrator_factory_get_instance(
    GauXCStatus* status,
    const GauXCIntegratorFactory integrator_factory,
    const GauXCFunctional func,
    const GauXCLoadBalancer lb
);

/**
 * @brief Get a shared XCIntegrator instance from an XCIntegratorFactory.
 * @param status Status object to capture any errors.
 * @param integrator_factory Handle to the XCIntegratorFactory.
 * @param func Handle to the XCFunctional.
 * @param lb Handle to the LoadBalancer.
 * @return Handle to the created XCIntegrator.
 */
extern GauXCIntegrator gauxc_integrator_factory_get_shared_instance(
    GauXCStatus* status,
    const GauXCIntegratorFactory integrator_factory,
    const GauXCFunctional func,
    const GauXCLoadBalancer lb
);

/**
 * @brief Integrate the density matrix to get the number of electrons.
 * @param status Status object to capture any errors.
 * @param integrator Handle to the XCIntegrator.
 * @param density_matrix Density matrix container.
 * @param den Pointer to store the number of electrons.
 */
extern void gauxc_integrator_integrate_den(
    GauXCStatus* status,
    const GauXCIntegrator integrator,
    const GauXCMatrix density_matrix,
    double* den
);

/**
 * @brief Evaluate the exchange-correlation energy for RKS.
 * @param status Status object to capture any errors.
 * @param integrator Handle to the XCIntegrator.
 * @param density_matrix Density matrix container for RKS.
 * @param exc Pointer to store the exchange-correlation energy.
 */
extern void gauxc_integrator_eval_exc_rks(
    GauXCStatus* status,
    const GauXCIntegrator integrator,
    const GauXCMatrix density_matrix,
    double* exc
);

/**
 * @brief Evaluate the exchange-correlation energy for UKS.
 * @param status Status object to capture any errors.
 * @param integrator Handle to the XCIntegrator.
 * @param density_matrix_s Density matrix container for total density.
 * @param density_matrix_z Density matrix container for spin density.
 * @param exc Pointer to store the exchange-correlation energy.
 */
extern void gauxc_integrator_eval_exc_uks(
    GauXCStatus* status,
    const GauXCIntegrator integrator,
    const GauXCMatrix density_matrix_s,
    const GauXCMatrix density_matrix_z,
    double* exc
);

/**
 * @brief Evaluate the exchange-correlation energy for GKS.
 * @param status Status object to capture any errors.
 * @param integrator Handle to the XCIntegrator.
 * @param density_matrix_s Density matrix container for total density.
 * @param density_matrix_z Density matrix container for spin z density.
 * @param density_matrix_x Density matrix container for spin x component.
 * @param density_matrix_y Density matrix container for spin y component.
 * @param exc Pointer to store the exchange-correlation energy.
 */
extern void gauxc_integrator_eval_exc_gks(
    GauXCStatus* status,
    const GauXCIntegrator integrator,
    const GauXCMatrix density_matrix_s,
    const GauXCMatrix density_matrix_z,
    const GauXCMatrix density_matrix_x,
    const GauXCMatrix density_matrix_y,
    double* exc
);

/**
 * @brief Evaluate the exchange-correlation energy and potential for RKS.
 * @param status Status object to capture any errors.
 * @param integrator Handle to the XCIntegrator.
 * @param density_matrix Density matrix container for RKS.
 * @param exc Pointer to store the exchange-correlation energy.
 * @param vxc_matrix Matrix container to store the exchange-correlation potential.
 */
extern void gauxc_integrator_eval_exc_vxc_rks(
    GauXCStatus* status,
    const GauXCIntegrator integrator,
    const GauXCMatrix density_matrix,
    double* exc,
    GauXCMatrix* vxc_matrix
);

/**
 * @brief Evaluate the exchange-correlation energy and potential for UKS.
 * @param status Status object to capture any errors.
 * @param integrator Handle to the XCIntegrator.
 * @param density_matrix_s Density matrix container for total density.
 * @param density_matrix_z Density matrix container for spin density.
 * @param exc Pointer to store the exchange-correlation energy.
 * @param vxc_matrix_s Matrix container to store the exchange-correlation potential for total density.
 * @param vxc_matrix_z Matrix container to store the exchange-correlation potential for spin density.
 */
extern void gauxc_integrator_eval_exc_vxc_uks(
    GauXCStatus* status,
    const GauXCIntegrator integrator,
    const GauXCMatrix density_matrix_s,
    const GauXCMatrix density_matrix_z,
    double* exc,
    GauXCMatrix* vxc_matrix_s,
    GauXCMatrix* vxc_matrix_z
);

/**
 * @brief Evaluate the exchange-correlation energy and potential for UKS.
 * @param status Status object to capture any errors.
 * @param integrator Handle to the XCIntegrator.
 * @param density_matrix_s Density matrix container for total density.
 * @param density_matrix_z Density matrix container for spin density.
 * @param model String specifying the OneDFT model to use.
 * @param exc Pointer to store the exchange-correlation energy.
 * @param vxc_matrix_s Matrix container to store the exchange-correlation potential for total density.
 * @param vxc_matrix_z Matrix container to store the exchange-correlation potential for spin density
 */
extern void gauxc_integrator_eval_exc_vxc_onedft_uks(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  const GauXCMatrix density_matrix_s,
  const GauXCMatrix density_matrix_z,
  const char* model,
  double* exc,
  GauXCMatrix* vxc_matrix_s,
  GauXCMatrix* vxc_matrix_z
);

/**
 * @brief Evaluate the exchange-correlation energy and potential for GKS.
 * @param status Status object to capture any errors.
 * @param integrator Handle to the XCIntegrator.
 * @param density_matrix_s Density matrix container for total density.
 * @param density_matrix_z Density matrix container for spin z density.
 * @param density_matrix_x Density matrix container for spin x component.
 * @param density_matrix_y Density matrix container for spin y component.
 * @param exc Pointer to store the exchange-correlation energy.
 * @param vxc_matrix_s Matrix container to store the exchange-correlation potential for total density.
 * @param vxc_matrix_z Matrix container to store the exchange-correlation potential for spin z density.
 * @param vxc_matrix_x Matrix container to store the exchange-correlation potential for spin x component.
 * @param vxc_matrix_y Matrix container to store the exchange-correlation potential for spin y component.
 */
extern void gauxc_integrator_eval_exc_vxc_gks(
    GauXCStatus* status,
    const GauXCIntegrator integrator,
    const GauXCMatrix density_matrix_s,
    const GauXCMatrix density_matrix_z,
    const GauXCMatrix density_matrix_x,
    const GauXCMatrix density_matrix_y,
    double* exc,
    GauXCMatrix* vxc_matrix_s,
    GauXCMatrix* vxc_matrix_z,
    GauXCMatrix* vxc_matrix_x,
    GauXCMatrix* vxc_matrix_y
);

#ifdef __cplusplus
} // extern "C"
} // namespace GauXC::C
#endif