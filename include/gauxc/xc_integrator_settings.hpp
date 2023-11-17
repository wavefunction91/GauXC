/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
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

struct IntegratorSettingsEXCVXC { virtual ~IntegratorSettingsEXCVXC() noexcept = default; };
struct IntegratorSettingsKS : public IntegratorSettingsEXCVXC {
  double gks_dtol = 1e-12;
};

}
