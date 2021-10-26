#ifndef _MY_BOYS_COMPUTATION
#define _MY_BOYS_COMPUTATION
#include <gauxc/util/constexpr_math.hpp>
#include <cmath>

namespace GauXC {

template <size_t M>
struct boys_asymp_const_coeff {

  static constexpr double value = 
    (constants::sqrt_pi<> / integral_pow_two<2*M+1>::value) *
    (integral_factorial<2*M>::value / integral_factorial<M>::value );

};

template <size_t M, typename T>
inline T boys_asymp( T x ) {
  const auto x_inv = 1./x;
  if constexpr (M == 0) {
    return constants::sqrt_pi_ov_2<> * std::sqrt( x_inv );
  } else {
    constexpr auto const_coeff = boys_asymp_const_coeff<M>::value;
    return const_coeff * std::sqrt( integral_pow<2*M+1>(x_inv) );
  }
  abort();
}

template <size_t M, typename T>
inline void boys_asymp( size_t npts, const T* X, T* FmX ) {
  #pragma unroll(NPTS_LOCAL)
  for( int i = 0; i < npts; ++i ) {
    FmX[i] = boys_asymp<M>(X[i]);
  }
}

template <size_t M>
double boys_function( double T );
template <size_t M>
void boys_function(int npts, const double* T, double* FmT );

}



double boys_asymp(int m, double T);
double boys_reference(int m, double T);
double boys_function(int m, double T);
void boys_function(int m, int npts, const double* T, double* FmT);
void boys_asymp(int npts, int m, const double* T, double* FmT);

#endif
