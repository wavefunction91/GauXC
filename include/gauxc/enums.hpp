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

#include <gauxc/enums.h>

namespace GauXC {

/**
 *  @brief GauXC specific enums for the specification of radial quadratures
 *
 *  Generally mapped to equivalent enums in IntegratorXX
 */
enum class RadialQuad {
  Becke = C::GauXC_RadialQuad_Becke,                         ///< Becke radial quadrature
  MuraKnowles = C::GauXC_RadialQuad_MuraKnowles,             ///< Mura-Knowles radial quadrature
  MurrayHandyLaming = C::GauXC_RadialQuad_MurrayHandyLaming, ///< Murray-Handy-Laming radial quadrature
  TreutlerAhlrichs = C::GauXC_RadialQuad_TreutlerAhlrichs    ///< Treutler-Ahlrichs radial quadrature
};

/**
 *  @brief Specifications of grid defaults for atomic integration
 *
 *  See https://gaussian.com/integral for specification
 */
enum class AtomicGridSizeDefault {
  FineGrid = C::GauXC_AtomicGridSizeDefault_FineGrid,            ///< Fine grid      (least accurate)
  UltraFineGrid = C::GauXC_AtomicGridSizeDefault_UltraFineGrid,  ///< Ultrafine grid (appropriate accuracy)
  SuperFineGrid = C::GauXC_AtomicGridSizeDefault_SuperFineGrid,  ///< Superfine grid (most accurate)
  GM3 = C::GauXC_AtomicGridSizeDefault_GM3,                      ///< Treutler-Ahlrichs GM3
  GM5 = C::GauXC_AtomicGridSizeDefault_GM5                       ///< Treutlet-Ahlrichs GM5
};

/**
 *  @brief Specifications of atomic partitioning scheme for the
 *  molecular integration
 */
enum class XCWeightAlg {
  NOTPARTITIONED = C::GauXC_XCWeightAlg_NOTPARTITIONED, ///< Not partitioned
  Becke = C::GauXC_XCWeightAlg_Becke,                   ///< The original Becke weighting scheme
  SSF = C::GauXC_XCWeightAlg_SSF,                       ///< The Stratmann-Scuseria-Frisch weighting scheme
  LKO = C::GauXC_XCWeightAlg_LKO                        ///< The Lauqua-Kuessman-Ochsenfeld weighting scheme
};

/**
 *  @brief Specification of the execution space for various operations
 */
enum class ExecutionSpace {
  Host = C::GauXC_ExecutionSpace_Host,    ///< Execute task on the host
  Device = C::GauXC_ExecutionSpace_Device ///< Execute task on the device (e.g. GPU)
};

/// Supported Algorithms / Integrands
enum class SupportedAlg {
  XC = C::GauXC_SupportedAlg_XC,
  DEN = C::GauXC_SupportedAlg_DEN,
  SNLINK = C::GauXC_SupportedAlg_SNLINK
};

/// High-level specification of pruning schemes for atomic quadratures
enum class PruningScheme {
  Unpruned = C::GauXC_PruningScheme_Unpruned, ///< Unpruned atomic quadrature
  Robust = C::GauXC_PruningScheme_Robust,     ///< The "Robust" scheme of Psi4
  Treutler = C::GauXC_PruningScheme_Treutler  ///< The Treutler-Ahlrichs scheme
};


} // namespace GauXC
