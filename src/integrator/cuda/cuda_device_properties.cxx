#include "cuda_device_properties.hpp"

namespace GauXC {
namespace cuda  {




uint32_t get_submat_cut_block(int32_t LDA) {
  return 512;
}

}
}
