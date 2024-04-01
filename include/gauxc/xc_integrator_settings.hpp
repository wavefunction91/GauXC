/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once

namespace GauXC {

struct IntegratorSettingsEXX { virtual ~IntegratorSettingsEXX() noexcept = default; };
struct IntegratorSettingsSNLinK : public IntegratorSettingsEXX {
  bool screen_ek = true;
  double energy_tol = 1e-10;
  double k_tol      = 1e-10;
};

struct IntegratorSettingsXC { virtual ~IntegratorSettingsXC() noexcept = default; };
struct IntegratorSettingsKS : public IntegratorSettingsXC {
  double gks_dtol = 1e-12;
};

}
