#pragma once

namespace GauXC {

/**
 *  @brief GauXC specific enums for the specification of radial quadratures
 *
 *  Generally mapped to equivalent enums in IntegratorXX
 */
enum class RadialQuad {
  MuraKnowles,      ///< Mura-Knowles radial quadrature
  MurrayHandyLaming ///< Murray-Handy-Laming radial quadrature
};

/**
 *  @brief Specifications of grid defaults for atomic integration
 *
 *  See https://gaussian.com/integral for specification
 */
enum class AtomicGridSizeDefault {
  FineGrid,       ///< Fine grid      (least accurate)
  UltraFineGrid,  ///< Ultrafine grid (appropriate accuracy)
  SuperFineGrid   ///< Superfine grid (most accurate)
};

/**
 *  @brief Specifications of atomic partitioning scheme for the
 *  molecular integration
 */
enum class XCWeightAlg {
  Becke, ///< The original Becke weighting scheme
  SSF,   ///< The Stratmann-Scuseria-Frisch weighting scheme
  LKO
};

/**
 *  @brief Specification of the execution space for various operations
 */
enum class ExecutionSpace {
  Host,  ///< Execute task on the host
  Device ///< Execute task on the device (e.g. GPU)
};


} // namespace GauXC
