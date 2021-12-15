#include "chebyshev_boys_computation.hpp"
#include <gauxc/util/constexpr_math.hpp>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <memory>
#include <vector>

#define MIN(a,b)			\
  ({ __typeof__ (a) _a = (a);	        \
  __typeof__ (b) _b = (b);		\
  _a < _b ? _a : _b; })

namespace GauXC {

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
  
  // create table - so this should be done on the host
  void cheby_coeff(int m, int ncheb, double a, double b, double* c) {
    const int n = ncheb+1;
    const double pi_ov_2n = M_PI / (2 * n);
    
    std::vector<double> f_table(n);
    for( int i = 0; i < n; ++i ) {
      double x = std::cos( (2.*(i+1)-1) * pi_ov_2n );
      x = 0.5 * ( a+b + (b-a)*x );
      f_table[i] = boys_reference(m, x);
    }

    c[0] = std::accumulate( f_table.begin(), f_table.end(),0. ) / n;
    for( int i = 1; i < n; ++i ) {
      double _val = 0.;
      for( int j = 0; j < n; ++j ) {
	_val += f_table[j] * std::cos( i * (2*(j+1)-1) * pi_ov_2n );
      }
      c[i] = 2.0 * _val / n;
    }
  }

  void cheby_to_monomial_coeff( int ncheb, double *coeff ) {
    const int n = ncheb+1;
    int64_t i_fact = 1;
    int64_t t_fact = 1;
    for(int i = 0; i < n; ++i) {
      if(i)     i_fact *= i;
      if(i > 1) t_fact *= 2;

      double _val = 0;
      if(!i) {
	int m1_fact = 1;
	for( int j = i; j < n; j += 2 ) {
	  _val += m1_fact * coeff[j];
	  m1_fact *= -1;
	}
      } else {
	int m1_term = 1;
	for( int j = i; j < n; j += 2 ) {
	  const int f_up = (i+j)/2 - 1;
	  const int f_lo = (j-i)/2;
	  int f_term = 1;
	  for( int k = f_lo+1; k <= f_up; ++k ) f_term *= k;
	  _val += t_fact * j * m1_term * double(f_term) / double(i_fact) * coeff[j];
	  m1_term *= -1;
	}

      }
      coeff[i] = _val;
    }
  }
  
  void generate_boys_table(int ncheb, int maxM, double maxT, int nseg, double* cheb_coeff_table, int ld) {
    const double deltaT = maxT / nseg;
    for( int m = 0; m <= maxM; ++m ) {
      double* coeff_m = cheb_coeff_table + m * ld * nseg; // table offset for current m
      for( int iseg = 0; iseg < nseg; ++iseg ) {
	double* coeff_seg = coeff_m + iseg * ld;

	const double a = iseg * deltaT;
	const double b = a + deltaT;

	cheby_coeff( m, ncheb, a, b, coeff_seg ); // Generate coeff in Chebyshev basis
	cheby_to_monomial_coeff( ncheb, coeff_seg );   // Convert to monomial basis
      }
    }
  }

  __device__ double *dev_boys_table_;
  double *host_boys_table_;
  
  void gauxc_boys_init() {
    double *tmp_ = (double*) malloc(DEFAULT_LD_TABLE * DEFAULT_NSEGMENT * (DEFAULT_MAX_M + 1) * sizeof(double));
    generate_boys_table(DEFAULT_NCHEB, DEFAULT_MAX_M, DEFAULT_MAX_T, DEFAULT_NSEGMENT, tmp_, DEFAULT_LD_TABLE);

    cudaMalloc((void**) &host_boys_table_, DEFAULT_LD_TABLE * DEFAULT_NSEGMENT * (DEFAULT_MAX_M + 1) * sizeof(double));
    cudaMemcpy(host_boys_table_, tmp_, DEFAULT_LD_TABLE * DEFAULT_NSEGMENT * (DEFAULT_MAX_M + 1) * sizeof(double), cudaMemcpyHostToDevice);
    
    cudaMemcpyToSymbol(dev_boys_table_, &host_boys_table_, sizeof(host_boys_table_));

    free(tmp_);
  }
  
  void gauxc_boys_finalize() {
    cudaFree(host_boys_table_);
  }

  // these functions are needed

  __device__ inline void monomial_expand(size_t npts, const double* coeff, const double *x, double a, double b, double* eval) {
    const int n = DEFAULT_NCHEB + 1;
    const double sum = a+b;
    const double diff = b-a;
    const double ratio = sum / diff;
    const double fact = 2. / diff;

    double xp[n]; xp[0] = 1.;

    for(size_t j = 0; j < npts; ++j) {
      double xt = fact * x[j] - ratio;

      for(int i = 1; i < n; ++i) xp[i] = xp[i-1] * xt;

      double _val = 0.;
      for(int i = 0; i < n; ++i) _val += xp[i] * coeff[i];

      eval[j] = _val;
    }
  }

  template <int M>
  __device__ inline double boys_asymp_element( double x ) {
    const auto x_inv = 1./x;

    if constexpr (M != 0) {
      constexpr double const_coeff = (constants::sqrt_pi<> / integral_pow_two<2*M+1>::value) * (integral_factorial<2*M>::value / integral_factorial<M>::value);
      return const_coeff * std::sqrt(integral_pow<2*M+1>(x_inv));
    }

    return constants::sqrt_pi_ov_2<> * std::sqrt( x_inv ); 
  }

  template <int M>
  __device__ inline void boys_asymp_elements( size_t npts, const double* X, double* FmX ) {
    for( size_t i = 0; i < npts; ++i ) {
      FmX[i] = boys_asymp_element<M>(X[i]);
    }
  }
  
  template <int M>
  __device__ double gauxc_boys_element(double T) {
    if constexpr (M != 0) {
      double* boys_m = (dev_boys_table_ + M * DEFAULT_LD_TABLE * DEFAULT_NSEGMENT);
      double deltaT = double(DEFAULT_MAX_T) / DEFAULT_NSEGMENT;

      if(T > DEFAULT_MAX_T) {
	return boys_asymp_element<M>(T);
      } else {
	double eval;
	  
	int iseg = std::floor(T/ deltaT);
	const double* boys_seg = boys_m + iseg * DEFAULT_LD_TABLE;

	const double a = iseg * deltaT;
	const double b = a + deltaT;
	monomial_expand(1, boys_seg, &T, a, b, &eval);

	return eval;
      }
    }

    return boys_asymp_element<M>(T);
  }
  
  template <int M>
  __device__ void gauxc_boys_elements(size_t npts, double* T, double* eval) {
    // Get min T
    auto min_t = std::min_element( T, T + npts );

    if( *min_t >  DEFAULT_MAX_T ) {
      boys_asymp_elements<M>(npts, T, eval);
    } else {
      if constexpr (M == 0) {
	  for(size_t i = 0; i < npts; ++i) {
	    const double sqrt_t = std::sqrt(T[i]);
	    const double inv_sqrt_t = 1./sqrt_t;
	    eval[i] = constants::sqrt_pi_ov_2<> * std::erf(sqrt_t) * inv_sqrt_t;
	  }
	} else {
	double* boys_m = dev_boys_table_ + M * DEFAULT_LD_TABLE * DEFAULT_NSEGMENT;
	double deltaT = double(DEFAULT_MAX_T) / DEFAULT_NSEGMENT;

	for( size_t i = 0; i < npts; ++i ) {
	  const double tval = T[i];
	  if( tval > DEFAULT_MAX_T ) {
	    eval[i] = boys_asymp_element<M>(tval); 
	  } else {
	    int iseg = std::floor(tval / deltaT);
	    const double* boys_seg = boys_m + iseg * DEFAULT_LD_TABLE;

	    const double a = iseg * deltaT;
	    const double b = a + deltaT;
	    monomial_expand(1, boys_seg, T+i, a, b, eval+i);
	  }
	}
      }
    }
  }

#define BOYS_FUNCTION_IMPLEMENTATION(M)					\
  template __device__ double gauxc_boys_element<M>(double);		\
  template __device__ void gauxc_boys_elements<M>(size_t, double*, double*);

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
  
  /*
  double boys_element(int m, double T) {
    switch(m) {
    case 0:
      return chebyshev_element<0>(T);
      break;
    case 1:
      return chebyshev_element<1>(T);
      break;
    case 2:
      return chebyshev_element<2>(T);
      break;
    case 3:
      return chebyshev_element<3>(T);
      break;
    case 4:
      return chebyshev_element<4>(T);
      break;
    case 5:
      return chebyshev_element<5>(T);
      break;
    case 6:
      return chebyshev_element<6>(T);
      break;
    case 7:
      return chebyshev_element<7>(T);
      break;
    case 8:
      return chebyshev_element<8>(T);
      break;
    case 9:
      return chebyshev_element<9>(T);
      break;
    case 10:
      return chebyshev_element<10>(T);
      break;
    case 11:
      return chebyshev_element<11>(T);
      break;
    case 12:
      return chebyshev_element<12>(T);
      break;
    case 13:
      return chebyshev_element<13>(T);
      break;
    case 14:
      return chebyshev_element<14>(T);
      break;
    case 15:
      return chebyshev_element<15>(T);
      break;
    default:
      return chebyshev_element<16>(T);
      break;
    }
  }

  void boys_elements(int m, size_t npts, double *T, double *eval) {
    switch(m) {
    case 0:
      chebyshev_elements<0>(npts, T, eval);
      break;
    case 1:
      chebyshev_elements<1>(npts, T, eval);
      break;
    case 2:
      chebyshev_elements<2>(npts, T, eval);
      break;
    case 3:
      chebyshev_elements<3>(npts, T, eval);
      break;
    case 4:
      chebyshev_elements<4>(npts, T, eval);
      break;
    case 5:
      chebyshev_elements<5>(npts, T, eval);
      break;
    case 6:
      chebyshev_elements<6>(npts, T, eval);
      break;
    case 7:
      chebyshev_elements<7>(npts, T, eval);
      break;
    case 8:
      chebyshev_elements<8>(npts, T, eval);
      break;
    case 9:
      chebyshev_elements<9>(npts, T, eval);
      break;
    case 10:
      chebyshev_elements<10>(npts, T, eval);
      break;
    case 11:
      chebyshev_elements<11>(npts, T, eval);
      break;
    case 12:
      chebyshev_elements<12>(npts, T, eval);
      break;
    case 13:
      chebyshev_elements<13>(npts, T, eval);
      break;
    case 14:
      chebyshev_elements<14>(npts, T, eval);
      break;
    case 15:
      chebyshev_elements<15>(npts, T, eval);
      break;
    default:
      chebyshev_elements<16>(npts, T, eval);
      break;
    }
  }
  */
}
