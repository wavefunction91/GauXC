#pragma once
#include <cmath>
#include <cstddef>

namespace GauXC {
namespace util  {

// R(r) = r^l * \sum_i c_i * exp(-a_i*r^2)
template <typename T>
T gau_rad_eval( int32_t l, int32_t nprim, const T* alpha, const T* coeff, T r ) {
  T tmp = 0.;
  const T r2 = r*r;
  for( auto i = 0; i < nprim; ++i ) {
    tmp += coeff[i] * std::exp( -alpha[i] * r2 );
  }
  return std::pow(r,l) * tmp;
}

template <typename T>
T gau_rad_cutoff( int32_t l, int32_t nprim, const T* alpha, const T* coeff, T tol ) {

  const double log_tol = -std::log(tol);
  // Initial guess
  double r = 0;
  for( auto i = 0; i < nprim; ++i ) {
    // Prim cutoff
    const double log_alpha = std::log(alpha[i]);
    const double prim_cutoff = 
      std::sqrt( (log_tol + log_alpha/2.)/alpha[i] );
    r = std::max( r, prim_cutoff );
  }

  std::vector<T> abs_coeff( coeff, coeff + nprim );
  for( auto& x : abs_coeff ) x = std::abs(x);
  double rad_eval = gau_rad_eval(l, nprim, alpha, abs_coeff.data(), r);

  const double step = 0.01;
  if( rad_eval > tol ) { 
    // Walk to the left
    while( rad_eval > tol ) {
      r = r + step;
      rad_eval = gau_rad_eval(l, nprim, alpha, abs_coeff.data(), r);
    }
  } else {
    // Walk to the right
    while( rad_eval < tol ) {
      r = r - step;
      rad_eval = gau_rad_eval(l, nprim, alpha, abs_coeff.data(), r);
    }
    // Correct for the extra step
    r = r + step;
    rad_eval = gau_rad_eval(l, nprim, alpha, abs_coeff.data(), r);
  }

  return r;

}

}
}
