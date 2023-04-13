/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
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
