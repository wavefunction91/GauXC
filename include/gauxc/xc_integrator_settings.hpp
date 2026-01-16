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
 *  @brief Base settings for exact exchange (EXX) integration.
 */
struct IntegratorSettingsEXX { virtual ~IntegratorSettingsEXX() noexcept = default; };

/**
 *  @brief Settings for seminumerical linear exchange (sn-LinK) integration.
 *
 *  Controls screening and tolerance parameters for sn-LinK calculations.
 */
struct IntegratorSettingsSNLinK : public IntegratorSettingsEXX {
  bool screen_ek = true;       ///< Whether to apply EK screening.
  double energy_tol = 1e-10;   ///< Energy convergence tolerance.
  double k_tol      = 1e-10;   ///< Exchange matrix element tolerance.
};

/**
 *  @brief Base settings for exchange-correlation (XC) integration.
 */
struct IntegratorSettingsXC { virtual ~IntegratorSettingsXC() noexcept = default; };

/**
 *  @brief Settings for Kohn-Sham DFT integration.
 *
 *  Controls tolerances for generalized Kohn-Sham (GKS) calculations.
 */
struct IntegratorSettingsKS : public IntegratorSettingsXC {
  double gks_dtol = 1e-12;  ///< Tolerance for GKS density evaluation.
};

/**
 *  @brief Settings for XC gradient evaluation.
 *
 *  Controls whether grid weight derivatives are included in the gradient
 *  computation or if only the Hellmann-Feynman contribution is used.
 */
struct IntegratorSettingsEXC_GRAD : public IntegratorSettingsKS {
  bool include_weight_derivatives = true; ///< Include grid weight contribution and translational invariance.
};

}
