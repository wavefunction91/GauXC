#include "boys_computation.h"
#include "chebyshev_boys_function.hpp"
#include <gauxc/util/constexpr_math.hpp>

#include <iostream>
#include <limits>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <array>
#include <iomanip>
#include <cassert>

double boys_asymp(int m, double T) {
  #if 1
  switch(m) {
    case 0:  return GauXC::boys_asymp<0> ( T );
    case 1:  return GauXC::boys_asymp<1> ( T );
    case 2:  return GauXC::boys_asymp<2> ( T );
    case 3:  return GauXC::boys_asymp<3> ( T );
    case 4:  return GauXC::boys_asymp<4> ( T );
    case 5:  return GauXC::boys_asymp<5> ( T );
    case 6:  return GauXC::boys_asymp<6> ( T );
    case 7:  return GauXC::boys_asymp<7> ( T );
    case 8:  return GauXC::boys_asymp<8> ( T );
    case 9:  return GauXC::boys_asymp<9> ( T );
    case 10: return GauXC::boys_asymp<10>( T );
    case 11: return GauXC::boys_asymp<11>( T );
    case 12: return GauXC::boys_asymp<12>( T );
    case 13: return GauXC::boys_asymp<13>( T );
    case 14: return GauXC::boys_asymp<14>( T );
    case 15: return GauXC::boys_asymp<15>( T );
    case 16: return GauXC::boys_asymp<16>( T );
  }
  abort();
  #else
  const auto one_ov_t = 1./T;
  const auto rsqrt_t  = std::sqrt(one_ov_t);
  double Fm = 0.88622692545275801365 * rsqrt_t;
  for (int i = 1; i <= m; ++i) {
    Fm = Fm * (i - 0.5) * one_ov_t;
  }
  return Fm;
  #endif
}

void boys_asymp(int npts, int m, const double* T, double* FmT) {

  switch(m) {
    case 0:  GauXC::boys_asymp<0> ( npts, T, FmT ); return;
    case 1:  GauXC::boys_asymp<1> ( npts, T, FmT ); return;
    case 2:  GauXC::boys_asymp<2> ( npts, T, FmT ); return;
    case 3:  GauXC::boys_asymp<3> ( npts, T, FmT ); return;
    case 4:  GauXC::boys_asymp<4> ( npts, T, FmT ); return;
    case 5:  GauXC::boys_asymp<5> ( npts, T, FmT ); return;
    case 6:  GauXC::boys_asymp<6> ( npts, T, FmT ); return;
    case 7:  GauXC::boys_asymp<7> ( npts, T, FmT ); return;
    case 8:  GauXC::boys_asymp<8> ( npts, T, FmT ); return;
    case 9:  GauXC::boys_asymp<9> ( npts, T, FmT ); return;
    case 10: GauXC::boys_asymp<10>( npts, T, FmT ); return;
    case 11: GauXC::boys_asymp<11>( npts, T, FmT ); return;
    case 12: GauXC::boys_asymp<12>( npts, T, FmT ); return;
    case 13: GauXC::boys_asymp<13>( npts, T, FmT ); return;
    case 14: GauXC::boys_asymp<14>( npts, T, FmT ); return;
    case 15: GauXC::boys_asymp<15>( npts, T, FmT ); return;
    case 16: GauXC::boys_asymp<16>( npts, T, FmT ); return;
  }

}

double boys_reference(int m, double T) {
  double denom = m + 0.5;
  double term  = std::exp(-T) / (2 * denom);
  double old_term = term;
  double sum = old_term;

  double eps = std::numeric_limits<double>::epsilon();
  double eps_10 = eps / 10;

  while( term > sum * eps_10 || old_term < term ) {
    denom = denom + 1;
    old_term = term;
    term = old_term * T / denom;
    sum = sum + term;
  }

  return sum;
}

double boys_function(int m, double T) {
  assert( GauXC::chebyshev_boys_instance );
  double val;
  GauXC::chebyshev_boys_instance->eval( 1, m, &T, &val );
  return val;
}

void boys_function(int m, int npts, const double* T, double* FmT) {
  assert( GauXC::chebyshev_boys_instance );
  GauXC::chebyshev_boys_instance->eval( npts, m, T, FmT );
}
