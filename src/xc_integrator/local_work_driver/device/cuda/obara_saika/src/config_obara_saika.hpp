#pragma once

#include <gauxc/util/constexpr_math.hpp>

#define NPTS_LOCAL 64

#define DEFAULT_NCHEB  7
#define DEFAULT_MAX_M  8
#define DEFAULT_MAX_T 30

#define DEFAULT_NSEGMENT ((DEFAULT_MAX_T * DEFAULT_NCHEB) / 2)
#define DEFAULT_LD_TABLE (DEFAULT_NCHEB + 1)
 
namespace XGPU {
  template <int M>
  __device__ inline void boys_element(double *T, double *T_inv_e, double *eval, double *boys_table) {
    if((*T) < DEFAULT_MAX_T) {
      if constexpr (M == 0) {
	const double sqrt_t = std::sqrt((*T));
	const double inv_sqrt_t = 1./sqrt_t;
	*(T_inv_e) = 0.0;
	*(eval) = GauXC::constants::sqrt_pi_ov_2<> * std::erf(sqrt_t) * inv_sqrt_t;
      } else {
	const double* boys_m = (boys_table + M * DEFAULT_LD_TABLE * DEFAULT_NSEGMENT);
	constexpr double deltaT = double(DEFAULT_MAX_T) / DEFAULT_NSEGMENT;
	constexpr double one_over_deltaT = 1 / deltaT;
	
	int iseg = std::floor((*T) * one_over_deltaT);
	const double* boys_seg = boys_m + iseg * DEFAULT_LD_TABLE;
	
	const double ratio = (2 * iseg + 1);
	const double fact  = 2.0 / deltaT;
	
	double xt = (*T) * fact - ratio;
	double _rec = 1.0;
	double _val = boys_seg[0];
	
	for(int i = 1; i < DEFAULT_NCHEB + 1; ++i) {
	  _rec = _rec * xt;
	  _val += _rec * boys_seg[i];
	}

	*(T_inv_e) = 0.5 * std::exp(-(*T));
	*(eval) = _val;
      }
    } else {
      const double t_inv = 1./(*T);
      double _val = GauXC::constants::sqrt_pi_ov_2<> * std::sqrt(t_inv);
    
      for(int i = 1; i < M + 1; ++i) {
	_val *= ((i - 0.5) * t_inv);
      }

      *(T_inv_e) = 0.0;
      *(eval) = _val;
    }
  }
}

#define SCALAR_TYPE double

#define SCALAR_LENGTH 1

#define SCALAR_SET1(x) (x)
#define SCALAR_ZERO() (0.0)

#define SCALAR_LOAD(x) *(x)
#define SCALAR_STORE(x, y) *(x) = y

#define SCALAR_ADD(x, y) (x + y)
#define SCALAR_SUB(x, y) (x - y)

#define SCALAR_MUL(x, y) (x * y)
#define SCALAR_FMA(x, y, z) (z + x * y)
#define SCALAR_FNMA(x, y, z) (z - x * y)

#define SCALAR_RECIPROCAL(x) (1.0 / (1.0 * x))

/*
  __device__  inline double monomial_expand(const double* coeff, const double x, double a, double b) {
  //const int n = DEFAULT_NCHEB + 1;
  const double sum = a+b;
  const double diff = b-a;
  const double ratio = sum / diff;
  const double fact = 2. / diff;

  //double xp[n]; xp[0] = 1.;
  double xp[DEFAULT_NCHEB + 1]; xp[0] = 1.;

  double xt = fact * x - ratio;

  //for(int i = 1; i < n; ++i) xp[i] = xp[i-1] * xt;
  for(int i = 1; i < DEFAULT_NCHEB + 1; ++i) xp[i] = xp[i-1] * xt;

  double _val = 0.;
  //for(int i = 0; i < n; ++i) _val += xp[i] * coeff[i];
  for(int i = 0; i < DEFAULT_NCHEB + 1; ++i) _val += xp[i] * coeff[i];

  return _val;
  }

  template <int M>
  __device__  inline double boys_asymp_element( double x ) {
  const auto x_inv = 1./x;

  if constexpr (M != 0) {
  constexpr double const_coeff = (constants::sqrt_pi<> / integral_pow_two<2*M+1>::value) * (integral_factorial<2*M>::value / integral_factorial<M>::value);
  return const_coeff * std::sqrt(integral_pow<2*M+1>(x_inv));
  }
    
  return constants::sqrt_pi_ov_2<> * std::sqrt( x_inv ); 
  }
  
  template <int M>
  __device__  inline double gauxc_boys_element(double *boys_table, double T) {

  if(T < DEFAULT_MAX_T) {
  if constexpr (M != 0) {
  const double* boys_m = (boys_table + M * DEFAULT_LD_TABLE * DEFAULT_NSEGMENT);
  constexpr double deltaT = double(DEFAULT_MAX_T) / DEFAULT_NSEGMENT;

  int iseg = std::floor(T/ deltaT);
  const double* boys_seg = boys_m + iseg * DEFAULT_LD_TABLE;

  const double a = iseg * deltaT;
  const double b = a + deltaT;
  return monomial_expand(boys_seg, T, a, b);
  }

  const double sqrt_t = std::sqrt(T);
  const double inv_sqrt_t = 1./sqrt_t;
  return constants::sqrt_pi_ov_2<> * std::erf(sqrt_t) * inv_sqrt_t;
  }

  return boys_asymp_element<M>(T);
  }
*/
