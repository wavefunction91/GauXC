#pragma once

namespace GauXC {

struct XCIntegratorState {
  bool load_balancer_populated     = false;
  bool modified_weights_are_stored = false;
};

}
