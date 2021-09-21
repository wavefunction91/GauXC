#pragma once

#include "ut_common.hpp"
#include <gauxc/molgrid.hpp>
#include <gauxc/basisset.hpp>
#include <gauxc/load_balancer.hpp>
#include <fstream>
#include <string>

using namespace GauXC;

#define MAX_NPTS_CHECK 67

struct ref_collocation_data {
  std::vector<int32_t>              mask;
  std::vector<std::array<double,3>> pts;
  std::vector<double>               eval;
  std::vector<double>               deval_x;
  std::vector<double>               deval_y;
  std::vector<double>               deval_z;

  template <typename Archive>
  void serialize( Archive& ar ) {
    ar( mask, pts, eval, deval_x, deval_y, deval_z );
  }

};

void check_collocation_transpose( int npts, int nbf, const double* ref_val, const double* comp_val, std::string msg = "" ) {

  // Check transpose
  for( int i = 0; i < nbf;  ++i )
  for( int j = 0; j < npts; ++j ) {
    INFO(msg << " IBF = " << i << " IPT = " << j);
    CHECK( ref_val[ i + j*nbf ] == Approx( comp_val[ i*npts + j ] ) );
  }

}

void check_collocation( int npts, int nbf, const double* ref_val, const double* comp_val ) {

  for( int i = 0; i < nbf;  ++i )
  for( int j = 0; j < npts; ++j ) {
    INFO("IBF = " << i << " IPT = " << j);
    CHECK( ref_val[ i + j*nbf ] == Approx( comp_val[ i + j*nbf ] ) );
  }

}
