#ifndef _MY_SIMD_INSTRUCTIONS
#define _MY_SIMD_INSTRUCTIONS

#include <gauxc/util/constexpr_math.hpp>

#define NPTS_LOCAL 64

#define SCALAR_TYPE double

#define SCALAR_LENGTH 1

#define SCALAR_SET1(x) x
#define SCALAR_ZERO() 0.0

#define SCALAR_LOAD(x) *(x)
#define SCALAR_STORE(x, y) *(x) = y

#define SCALAR_ADD(x, y) x + y
#define SCALAR_SUB(x, y) x - y

#define SCALAR_MUL(x, y) x * y
#define SCALAR_FMA(x, y, z) z + x * y 
#define SCALAR_FNMA(x, y, z) z - x * y

#define SCALAR_RECIPROCAL(x) 1.0 / (1.0 * x)

// these functions are needed

#define DEFAULT_NCHEB  7
#define DEFAULT_MAX_M 16
#define DEFAULT_MAX_T 30

#define DEFAULT_NSEGMENT ((DEFAULT_MAX_T * DEFAULT_NCHEB) / 2)
#define DEFAULT_LD_TABLE (DEFAULT_NCHEB + 1)
 
namespace GauXC {

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
    if constexpr (M != 0) {
	if (T < DEFAULT_MAX_T) {
	  double* boys_m = (boys_table + M * DEFAULT_LD_TABLE * DEFAULT_NSEGMENT);
	  double deltaT = double(DEFAULT_MAX_T) / DEFAULT_NSEGMENT;
	  
	  int iseg = std::floor(T/ deltaT);
	  const double* boys_seg = (boys_m + iseg * DEFAULT_LD_TABLE);
	  
	  const double a = iseg * deltaT;
	  const double b = a + deltaT;
	  
	  return monomial_expand(boys_seg, T, a, b);
	}
      }

    return boys_asymp_element<M>(T);
  }

#define BOYS_FUNCTION_IMPLEMENTATION(M)					\
  template __device__ double gauxc_boys_element<M>(double *, double);	\

  BOYS_FUNCTION_IMPLEMENTATION( 0);
  BOYS_FUNCTION_IMPLEMENTATION( 1);
  BOYS_FUNCTION_IMPLEMENTATION( 2);
  BOYS_FUNCTION_IMPLEMENTATION( 3);
  BOYS_FUNCTION_IMPLEMENTATION( 4);
  BOYS_FUNCTION_IMPLEMENTATION( 5);
  BOYS_FUNCTION_IMPLEMENTATION( 6);
  BOYS_FUNCTION_IMPLEMENTATION( 7);
  BOYS_FUNCTION_IMPLEMENTATION( 8);
  BOYS_FUNCTION_IMPLEMENTATION( 9);
  BOYS_FUNCTION_IMPLEMENTATION(10);
  BOYS_FUNCTION_IMPLEMENTATION(11);
  BOYS_FUNCTION_IMPLEMENTATION(12);
  BOYS_FUNCTION_IMPLEMENTATION(13);
  BOYS_FUNCTION_IMPLEMENTATION(14);
  BOYS_FUNCTION_IMPLEMENTATION(15);
}

#endif
