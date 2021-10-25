#include "boys_computation.h"
#include "chebyshev_boys_function.hpp"

#include <iostream>
#include <limits>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <array>
#include <iomanip>
#include <cassert>

#define PI           3.14159265358979323846
#define SQRT_PI      1.77245385090551602729
#define SQRT_PI_OV_2 0.88622692545275801364

#define MIN(a,b)			\
  ({ __typeof__ (a) _a = (a);	        \
  __typeof__ (b) _b = (b);		\
  _a < _b ? _a : _b; })


int64_t difact( int64_t i ) {
  int64_t v = 1;
  for( int k = 0; k < (i/2); ++k ) v *= i - 2 * k;
  return v;
}

template <size_t N>
struct pow_two : std::integral_constant< size_t, (1ul << N) > { };

template <size_t N>
struct integral_factorial;

template <>
struct integral_factorial<0ul> : std::integral_constant< size_t, 1ul > { };
template <size_t N>
struct integral_factorial : std::integral_constant< size_t, N * integral_factorial<N-1>::value > { };


template <size_t N, typename T>
constexpr T integral_pow(T x) {
  if constexpr (N == 0) { return T(1.); }
  else if constexpr (N == 1 ) { return x; }
  else { return x * integral_pow<N-1>(x); }
}

template <size_t M>
struct boys_asymp_const_coeff {

  static constexpr double value = (SQRT_PI / pow_two<2*M+1>::value) *
                                  (integral_factorial<2*M>::value / integral_factorial<M>::value );

};

template <size_t M, typename T>
T boys_asymp_eval( T x ) {
  constexpr auto const_coeff = boys_asymp_const_coeff<M>::value;
  const auto x_inv = 1./x;
  return const_coeff * std::sqrt( integral_pow<2*M+1>(x_inv) );
}

double boys_asymp(int m, double T) {
  #if 1
  switch(m) {
    case 0:  return SQRT_PI_OV_2 * sqrt( 1. / T );
    case 1:  return boys_asymp_eval<1> ( T );
    case 2:  return boys_asymp_eval<2> ( T );
    case 3:  return boys_asymp_eval<3> ( T );
    case 4:  return boys_asymp_eval<4> ( T );
    case 5:  return boys_asymp_eval<5> ( T );
    case 6:  return boys_asymp_eval<6> ( T );
    case 7:  return boys_asymp_eval<7> ( T );
    case 8:  return boys_asymp_eval<8> ( T );
    case 9:  return boys_asymp_eval<9> ( T );
    case 10: return boys_asymp_eval<10>( T );
    case 11: return boys_asymp_eval<11>( T );
    case 12: return boys_asymp_eval<12>( T );
    case 13: return boys_asymp_eval<13>( T );
    case 14: return boys_asymp_eval<14>( T );
    case 15: return boys_asymp_eval<15>( T );
    case 16: return boys_asymp_eval<16>( T );
  }
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

  for(int i_st = 0; i_st < npts; i_st += NPTS_LOCAL) {
    int ndo = MIN( NPTS_LOCAL, npts - i_st );
    #pragma unroll NPTS_LOCAL
    for( int i = 0; i < ndo; ++i ) {
      switch(m) {
        case 0:  FmT[i] = SQRT_PI_OV_2 * sqrt( 1. / T[i] );
        case 1:  FmT[i] = boys_asymp_eval<1> ( T[i] );
        case 2:  FmT[i] = boys_asymp_eval<2> ( T[i] );
        case 3:  FmT[i] = boys_asymp_eval<3> ( T[i] );
        case 4:  FmT[i] = boys_asymp_eval<4> ( T[i] );
        case 5:  FmT[i] = boys_asymp_eval<5> ( T[i] );
        case 6:  FmT[i] = boys_asymp_eval<6> ( T[i] );
        case 7:  FmT[i] = boys_asymp_eval<7> ( T[i] );
        case 8:  FmT[i] = boys_asymp_eval<8> ( T[i] );
        case 9:  FmT[i] = boys_asymp_eval<9> ( T[i] );
        case 10: FmT[i] = boys_asymp_eval<10>( T[i] );
        case 11: FmT[i] = boys_asymp_eval<11>( T[i] );
        case 12: FmT[i] = boys_asymp_eval<12>( T[i] );
        case 13: FmT[i] = boys_asymp_eval<13>( T[i] );
        case 14: FmT[i] = boys_asymp_eval<14>( T[i] );
        case 15: FmT[i] = boys_asymp_eval<15>( T[i] );
        case 16: FmT[i] = boys_asymp_eval<16>( T[i] );
      }
    }
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
