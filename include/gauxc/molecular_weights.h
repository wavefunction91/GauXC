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
#include <gauxc/enums.h>
#include <gauxc/load_balancer.h>

#ifdef __cplusplus
namespace GauXC::C {
extern "C" {
#endif

/**
 * @brief GauXC C API MolecularWeightsSettings representation.
 */
typedef struct GauXCMolecularWeightsSettings {
  enum GauXC_XCWeightAlg weight_alg; ///< Weight partitioning scheme
  bool becke_size_adjustment;        ///< Whether to use Becke size adjustments
} GauXCMolecularWeightsSettings;

/**
 * @brief GauXC C API MolecularWeights handle.
 */
typedef struct GauXCMolecularWeights {
  GauXCHeader hdr; ///< Header for internal use.
  void* ptr; ///< Pointer to the MolecularWeights instance.
  bool owned; ///< Whether this instance owns the MolecularWeights.
} GauXCMolecularWeights;

/**
 * @brief Delete a MolecularWeights instance.
 * @param status Status object to capture any errors.
 * @param mw Handle to the MolecularWeights to delete.
 */
extern void gauxc_molecular_weights_delete(
    GauXCStatus* status,
    GauXCMolecularWeights* mw
);

/**
 * @brief Apply molecular weights to a LoadBalancer's tasks.
 * @param status Status object to capture any errors.
 * @param mw Handle to the MolecularWeights.
 * @param lb Handle to the LoadBalancer.
 */
extern void gauxc_molecular_weights_modify_weights(
    GauXCStatus* status,
    const GauXCMolecularWeights mw,
    const GauXCLoadBalancer lb
);

/**
 * @brief GauXC C API MolecularWeightsFactory handle.
 */
typedef struct GauXCMolecularWeightsFactory {
  GauXCHeader hdr; ///< Header for internal use.
  void* ptr; ///< Pointer to the MolecularWeightsFactory instance.
} GauXCMolecularWeightsFactory;

/**
 * @brief Create a new MolecularWeightsFactory instance.
 * @param status Status object to capture any errors.
 * @param ex Execution space.
 * @param local_work_kernel_name Name of the LocalWorkDriver kernel to use.
 * @param settings Settings for the MolecularWeights calculation.
 * @return Handle to the created MolecularWeightsFactory.
 */
extern GauXCMolecularWeightsFactory gauxc_molecular_weights_factory_new(
    GauXCStatus* status,
    enum GauXC_ExecutionSpace ex,
    const char* local_work_kernel_name,
    GauXCMolecularWeightsSettings settings
);

/**
 * @brief Delete a MolecularWeightsFactory instance.
 * @param status Status object to capture any errors.
 * @param factory Handle to the MolecularWeightsFactory to delete.
 */
extern void gauxc_molecular_weights_factory_delete(
    GauXCStatus* status,
    GauXCMolecularWeightsFactory* factory
);

/**
 * @brief Get MolecularWeights instance from a MolecularWeightsFactory.
 * @param status Status object to capture any errors.
 * @param factory Handle to the MolecularWeightsFactory.
 * @return Handle to the created MolecularWeights.
 */
extern GauXCMolecularWeights gauxc_molecular_weights_factory_get_instance(
    GauXCStatus* status,
    const GauXCMolecularWeightsFactory factory
);

/**
 * @brief Get shared MolecularWeights instance from a MolecularWeightsFactory.
 * @param status Status object to capture any errors.
 * @param factory Handle to the MolecularWeightsFactory.
 * @return Handle to the created MolecularWeights.
 */
extern GauXCMolecularWeights gauxc_molecular_weights_factory_get_shared_instance(
    GauXCStatus* status,
    const GauXCMolecularWeightsFactory factory
);



#ifdef __cplusplus
} // extern "C"
} // namespace GauXC::C
#endif