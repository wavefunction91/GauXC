#pragma once
#include <vector>
#include <cstdint>

#include <gauxc/gauxc_config.hpp>

namespace GauXC {

template <typename F>
struct XCHostData {

  std::vector<F> eps;
  std::vector<F> gamma;
  std::vector<F> vrho;
  std::vector<F> vgamma;
 
  std::vector<F> zmat;
  std::vector<F> gmat;
  std::vector<F> nbe_scr;
  std::vector<F> den_scr;
  std::vector<F> basis_eval;
   
  inline XCHostData() {}

};

}
