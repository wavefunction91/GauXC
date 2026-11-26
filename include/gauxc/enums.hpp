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

namespace GauXC {

/**
 *  @brief GauXC specific enums for the specification of radial quadratures
 *
 *  Generally mapped to equivalent enums in IntegratorXX
 */
enum class RadialQuad {
  Becke,              ///< Becke radial quadrature
  MuraKnowles,        ///< Mura-Knowles radial quadrature
  MurrayHandyLaming,  ///< Murray-Handy-Laming radial quadrature
  TreutlerAhlrichs    ///< Treutler-Ahlrichs radial quadrature
};

/**
 *  @brief Specifications of grid defaults for atomic integration
 *
 *  See https://gaussian.com/integral for specification
 */
enum class AtomicGridSizeDefault {
  FineGrid,       ///< Fine grid      (least accurate)
  UltraFineGrid,  ///< Ultrafine grid (appropriate accuracy)
  SuperFineGrid,  ///< Superfine grid (most accurate)
  GM3,            ///< Treutler-Ahlrichs GM3
  GM5,            ///< Treutler-Ahlrichs GM5
  PySCF0,         ///< PySCF default level 0
  PySCF1,         ///< PySCF default level 1
  PySCF2,         ///< PySCF default level 2 (angular points ~ fine grid)
  PySCF3,         ///< PySCF default level 3
  PySCF4,         ///< PySCF default level 4 (radial points ~ fine grid, angular points ~ ultrafine grid)
  PySCF5,         ///< PySCF default level 5
  PySCF6,         ///< PySCF default level 6 (radial points ~ ultrafine grid, angular points ~ superfine grid)
  PySCF7,         ///< PySCF default level 7
  PySCF8,         ///< PySCF default level 8
  PySCF9          ///< PySCF default level 9 (radial points ~ superfine grid)
};

/**
 *  @brief Specifications of atomic partitioning scheme for the
 *  molecular integration
 */
enum class XCWeightAlg {
  NOTPARTITIONED,  ///< Not partitioned
  Becke,           ///< The original Becke weighting scheme
  SSF,             ///< The Stratmann-Scuseria-Frisch weighting scheme
  LKO              ///< The Lauqua-Kuessman-Ochsenfeld weighting scheme
};

/**
 *  @brief Specification of the execution space for various operations
 */
enum class ExecutionSpace {
  Host,   ///< Execute task on the host
  Device  ///< Execute task on the device (e.g. GPU)
};

/// Supported Algorithms / Integrands
enum class SupportedAlg {
  XC,
  DEN,
  SNLINK,
};

/// High-level specification of pruning schemes for atomic quadratures
enum class PruningScheme {
  Unpruned,       ///< Unpruned atomic quadrature
  Robust,         ///< The "Robust" scheme of Psi4
  Treutler,       ///< The Treutler-Ahlrichs scheme
  PySCF_Treutler, ///< The Treutler scheme of PySCF
  PySCF_SG1,      ///< The SG1 scheme of PySCF
  PySCF_NWChem,   ///< The NWChem scheme of PySCF
  PySCF_SGX       ///< The SGX scheme of PySCF
};


} // namespace GauXC
