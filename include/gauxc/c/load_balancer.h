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
#include <cstddef>
#include <cstdbool>
#else
#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#endif
#include <gauxc/c/types.h>
#include <gauxc/c/status.h>
#include <gauxc/c/enums.h>
#include <gauxc/c/runtime_environment.h>
#include <gauxc/c/molecule.h>
#include <gauxc/c/molgrid.h>
#include <gauxc/c/basisset.h>

#ifdef __cplusplus
extern "C" {
namespace GauXC::C {
#endif

/**
 * @brief GauXC C API LoadBalancer handle.
 */
typedef struct GauXCLoadBalancer {
  GauXCHeader hdr; ///< Header for internal use.
  void* ptr; ///< Pointer to the LoadBalancer instance.
} GauXCLoadBalancer;

/**
 * @brief Delete a LoadBalancer instance.
 * @param status Status object to capture any errors.
 * @param lb Handle to the LoadBalancer to delete.
 */
extern void gauxc_load_balancer_delete(
    GauXCStatus* status,
    GauXCLoadBalancer* lb
);

/**
 * @brief GauXC C API LoadBalancerFactory handle.
 */
typedef struct GauXCLoadBalancerFactory {
  GauXCHeader hdr; ///< Header for internal use.
  void* ptr; ///< Pointer to the LoadBalancerFactory instance.
} GauXCLoadBalancerFactory;

/**
 * @brief Create a new LoadBalancerFactory instance.
 * @param status Status object to capture any errors.
 * @param ex Execution space.
 * @param kernel_name Name of the load balancing kernel to use.
 * @return Handle to the created LoadBalancerFactory.
 */
extern GauXCLoadBalancerFactory gauxc_load_balancer_factory_new(
    GauXCStatus* status,
    enum GauXC_ExecutionSpace ex,
    const char* kernel_name
);

/**
 * @brief Delete a LoadBalancerFactory instance.
 * @param status Status object to capture any errors.
 * @param factory Handle to the LoadBalancerFactory to delete.
 */
extern void gauxc_load_balancer_factory_delete(
    GauXCStatus* status,
    GauXCLoadBalancerFactory* factory
);

/**
 * @brief Create a new LoadBalancer instance from a LoadBalancerFactory.
 * @param status Status object to capture any errors.
 * @param factory Handle to the LoadBalancerFactory.
 * @param env Handle to the RuntimeEnvironment.
 * @param mol Handle to the Molecule.
 * @param mg Handle to the MolGrid.
 * @param basis Handle to the BasisSet.
 * @return Handle to the created LoadBalancer.
 */
extern GauXCLoadBalancer gauxc_load_balancer_factory_get_instance(
    GauXCStatus* status,
    const GauXCLoadBalancerFactory factory,
    const GauXCRuntimeEnvironment env,
    const GauXCMolecule mol,
    const GauXCMolGrid mg,
    const GauXCBasisSet basis
);

#ifdef __cplusplus
} // namespace GauXC::C
} // extern "C"
#endif