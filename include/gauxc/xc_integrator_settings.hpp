#pragma once

namespace GauXC {

struct IntegratorSettingsEXX { virtual ~IntegratorSettingsEXX() noexcept = default; };
struct IntegratorSettingsSNLinK : public IntegratorSettingsEXX {
  bool screen_ek = true;
  double energy_tol = 1e-10;
  double k_tol      = 1e-10;
};

}
