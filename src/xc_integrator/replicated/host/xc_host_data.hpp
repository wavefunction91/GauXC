/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
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
  std::vector<F> tau;
  std::vector<F> lapl;
  std::vector<F> vrho;
  std::vector<F> vgamma;
  std::vector<F> vtau;
  std::vector<F> vlapl;
 
  std::vector<F> zmat;
  std::vector<F> gmat;
  std::vector<F> nbe_scr;
  std::vector<F> den_scr;
  std::vector<F> basis_eval;

  std::vector<F> epc;
  std::vector<F> protonic_vrho;
 
  std::vector<F> protonic_zmat;
  std::vector<F> protonic_gmat;
  std::vector<F> protonic_den_scr;
  std::vector<F> protonic_basis_eval;
   
  inline XCHostData() {}

};

}
