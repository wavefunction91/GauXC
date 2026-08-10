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
extern "C" {
namespace GauXC::C {
#endif

/**
 *  @brief GauXC specific enums for the specification of radial quadratures
 *
 *  Generally mapped to equivalent enums in IntegratorXX
 */
enum GauXC_RadialQuad {
  GauXC_RadialQuad_Becke,             ///< Becke radial quadrature
  GauXC_RadialQuad_MuraKnowles,       ///< Mura-Knowles radial quadrature
  GauXC_RadialQuad_MurrayHandyLaming, ///< Murray-Handy-Laming radial quadrature
  GauXC_RadialQuad_TreutlerAhlrichs   ///< Treutler-Ahlrichs radial quadrature
};

/**
 *  @brief Specifications of grid defaults for atomic integration
 *
 *  See https://gaussian.com/integral for specification
 */
enum GauXC_AtomicGridSizeDefault {
  GauXC_AtomicGridSizeDefault_FineGrid,       ///< Fine grid      (least accurate)
  GauXC_AtomicGridSizeDefault_UltraFineGrid,  ///< Ultrafine grid (appropriate accuracy)
  GauXC_AtomicGridSizeDefault_SuperFineGrid,  ///< Superfine grid (most accurate)
  GauXC_AtomicGridSizeDefault_GM3,            ///< Treutler-Ahlrichs GM3
  GauXC_AtomicGridSizeDefault_GM5,            ///< Treutler-Ahlrichs GM5
  GauXC_AtomicGridSizeDefault_PySCF0,         ///< PySCF default level 0
  GauXC_AtomicGridSizeDefault_PySCF1,         ///< PySCF default level 1
  GauXC_AtomicGridSizeDefault_PySCF2,         ///< PySCF default level 2 (angular points ~ fine grid)
  GauXC_AtomicGridSizeDefault_PySCF3,         ///< PySCF default level 3
  GauXC_AtomicGridSizeDefault_PySCF4,         ///< PySCF default level 4 (radial points ~ fine grid, angular points ~ ultrafine grid)
  GauXC_AtomicGridSizeDefault_PySCF5,         ///< PySCF default level 5
  GauXC_AtomicGridSizeDefault_PySCF6,         ///< PySCF default level 6 (radial points ~ ultrafine grid, angular points ~ superfine grid)
  GauXC_AtomicGridSizeDefault_PySCF7,         ///< PySCF default level 7
  GauXC_AtomicGridSizeDefault_PySCF8,         ///< PySCF default level 8
  GauXC_AtomicGridSizeDefault_PySCF9          ///< PySCF default level 9 (radial points ~ superfine grid)
};

/**
 *  @brief Specifications of atomic partitioning scheme for the
 *  molecular integration
 */
enum GauXC_XCWeightAlg {
  GauXC_XCWeightAlg_NOTPARTITIONED, ///< Not partitioned
  GauXC_XCWeightAlg_Becke, ///< The original Becke weighting scheme
  GauXC_XCWeightAlg_SSF,   ///< The Stratmann-Scuseria-Frisch weighting scheme
  GauXC_XCWeightAlg_LKO    ///< The Lauqua-Kuessman-Ochsenfeld weighting scheme
};

/**
 *  @brief Specification of the execution space for various operations
 */
enum GauXC_ExecutionSpace {
  GauXC_ExecutionSpace_Host,  ///< Execute task on the host
  GauXC_ExecutionSpace_Device ///< Execute task on the device (e.g. GPU)
};

/// Supported Algorithms / Integrands
enum GauXC_SupportedAlg {
  GauXC_SupportedAlg_XC,     ///< Exchange-correlation integration
  GauXC_SupportedAlg_DEN,    ///< Density integration
  GauXC_SupportedAlg_SNLINK  ///< Seminumerical Coulomb/exchange (snLinK)
};

/// High-level specification of pruning schemes for atomic quadratures
enum GauXC_PruningScheme {
  GauXC_PruningScheme_Unpruned,       ///< Unpruned atomic quadrature
  GauXC_PruningScheme_Robust,         ///< The "Robust" scheme of Psi4
  GauXC_PruningScheme_Treutler,       ///< The Treutler-Ahlrichs scheme
  GauXC_PruningScheme_PySCF_Treutler, ///< The Treutler scheme of PySCF
  GauXC_PruningScheme_PySCF_SG1,      ///< The SG1 scheme of PySCF
  GauXC_PruningScheme_PySCF_NWChem,   ///< The NWChem scheme of PySCF
  GauXC_PruningScheme_PySCF_SGX       ///< The SGX scheme of PySCF
};


#ifdef __cplusplus
} // namespace GauXC::C
} // extern "C"
#endif