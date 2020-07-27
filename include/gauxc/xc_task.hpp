#pragma once

#include <array>
#include <vector>
#include <cstdint>
#include "gauxc_config.hpp"

namespace GauXC {

struct XCTask {

  int32_t                             iParent;
  std::vector< std::array<double,3> > points;
  std::vector< double  >              weights;
  std::vector< int32_t >              shell_list;
  int32_t                             nbe;

  double                              dist_nearest;

  void merge_with( const XCTask& );


#ifdef GAUXC_ENABLE_CEREAL
  template <typename Archive>
  void serialize( Archive& ar ) {
    ar( iParent, nbe, dist_nearest, shell_list, points, weights );  
  }
#endif

};


}
