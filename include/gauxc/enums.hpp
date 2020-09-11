#pragma once

namespace GauXC {

enum class RadialQuad {
  MuraKnowles,
  MurrayHandyLaming
};

enum class AtomicGridSizeDefault {
  FineGrid,
  UltraFineGrid,
  SuperFineGrid
};

enum class XCWeightAlg {
  Becke,
  SSF
};

enum class ExecutionSpace {
  Host,
  Device
};


}
