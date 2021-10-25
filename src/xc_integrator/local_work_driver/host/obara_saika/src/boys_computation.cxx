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
struct pow_two;

#if 0
template<>
struct pow_two<0ul> : std::integral_constant<size_t, 1> { };

template <size_t N>
struct pow_two : std::integral_constant< size_t, 2ul * pow_two<N-1>::value > { };
#else
template <size_t N>
struct pow_two : std::integral_constant< size_t, (1ul << N) > { };
#endif

double my_pow_two( size_t m ) {
  switch(m) {
    case 0:  return pow_two<0 >::value;
    case 1:  return pow_two<1 >::value;
    case 2:  return pow_two<2 >::value;
    case 3:  return pow_two<3 >::value;
    case 4:  return pow_two<4 >::value;
    case 5:  return pow_two<5 >::value;
    case 6:  return pow_two<6 >::value;
    case 7:  return pow_two<7 >::value;
    case 8:  return pow_two<8 >::value;
    case 9:  return pow_two<9 >::value;
    case 10: return pow_two<10>::value;
    case 11: return pow_two<11>::value;
    case 12: return pow_two<12>::value;
    case 13: return pow_two<13>::value;
    case 14: return pow_two<14>::value;
    case 15: return pow_two<15>::value;
    default: return 0;
  }
}

template <size_t N>
struct integral_factorial;

template <>
struct integral_factorial<0ul> : std::integral_constant< size_t, 1ul > { };
template <size_t N>
struct integral_factorial : std::integral_constant< size_t, N * integral_factorial<N-1>::value > { };

#if 0
int64_t ifact( int64_t i ) {
  if((i == 0) || (i == 1)) return 1;
  int64_t v = 1;
  for( int k = 1; k <= i; ++k ) v *= k;
  return v;
}
#else
double ifact( size_t i ) {
  switch(i) {
    case 0:  return integral_factorial<0 >::value;
    case 1:  return integral_factorial<1 >::value;
    case 2:  return integral_factorial<2 >::value;
    case 3:  return integral_factorial<3 >::value;
    case 4:  return integral_factorial<4 >::value;
    case 5:  return integral_factorial<5 >::value;
    case 6:  return integral_factorial<6 >::value;
    case 7:  return integral_factorial<7 >::value;
    case 8:  return integral_factorial<8 >::value;
    case 9:  return integral_factorial<9 >::value;
    case 10: return integral_factorial<10>::value;
    case 11: return integral_factorial<11>::value;
    case 12: return integral_factorial<12>::value;
    case 13: return integral_factorial<13>::value;
    case 14: return integral_factorial<14>::value;
    case 15: return integral_factorial<15>::value;
    default: return 0;
  }
}
#endif

template <size_t N, typename T>
constexpr T integral_pow(T x) {
  if constexpr (N == 0) { return T(1.); }
  else if constexpr (N == 1 ) { return x; }
  else { return x * integral_pow<N-1>(x); }
}

template <typename T>
T ipow( T x, size_t n ) {
  switch(n) {
    case 0:  return integral_pow<0 >(x);
    case 1:  return integral_pow<1 >(x);
    case 2:  return integral_pow<2 >(x);
    case 3:  return integral_pow<3 >(x);
    case 4:  return integral_pow<4 >(x);
    case 5:  return integral_pow<5 >(x);
    case 6:  return integral_pow<6 >(x);
    case 7:  return integral_pow<7 >(x);
    case 8:  return integral_pow<8 >(x);
    case 9:  return integral_pow<9 >(x);
    case 10: return integral_pow<10>(x);
    case 11: return integral_pow<11>(x);
    case 12: return integral_pow<12>(x);
    case 13: return integral_pow<13>(x);
    case 14: return integral_pow<14>(x);
    case 15: return integral_pow<15>(x);
    default: return 0;
  }
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
  #if 0
  return difact(2 * m - 1) / 
         pow(2.,m + 1) *
         sqrt(PI / pow(T, 2 * m + 1));
  #else
  #if 0
  if( m == 0 ) {
    return SQRT_PI * 0.5 * sqrt( 1. / T );
  } else {
    return SQRT_PI * ifact(2*m) / my_pow_two(2*m+1) / ifact(m) / sqrt(ipow(T, 2 * m + 1));
  }
  #else
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
  #endif
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
#if 0
  if(T < 1e-10) {
    return 1. / (2*m + 1);
  } else if(T < 117) {
    return boys_reference(m, T);
  } else {
    return boys_asymp(m, T);
  }
#else
  assert( GauXC::chebyshev_boys_instance );
  double val;
  GauXC::chebyshev_boys_instance->eval( 1, m, &T, &val );
  return val;
#endif

}

void boys_function(int m, int npts, const double* T, double* FmT) {
  assert( GauXC::chebyshev_boys_instance );
  GauXC::chebyshev_boys_instance->eval( npts, m, T, FmT );
}
