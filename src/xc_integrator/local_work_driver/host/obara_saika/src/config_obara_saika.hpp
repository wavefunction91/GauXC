#pragma once

#include <gauxc/util/constexpr_math.hpp>

//#define X86_AVX
#define NPTS_LOCAL 64

#define DEFAULT_NCHEB  7
#define DEFAULT_MAX_M  8
#define DEFAULT_MAX_T 30

#define DEFAULT_NSEGMENT ((DEFAULT_MAX_T * DEFAULT_NCHEB) / 2)
#define DEFAULT_LD_TABLE (DEFAULT_NCHEB + 1)

namespace XCPU {
  template <int M>
  inline void boys_element(double *T, double *T_inv_e, double *eval, double *boys_table) {
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

  template <int M>
  inline void boys_elements(size_t npts, double* T, double *T_inv_e, double* eval, double *boys_table) {    
    for(size_t i = 0; i < npts; ++i) {
      if(T[i] < DEFAULT_MAX_T) {
	if constexpr (M == 0) {
	  const double sqrt_t = std::sqrt(T[i]);
	  const double inv_sqrt_t = 1./sqrt_t;
	  
	  T_inv_e[i] = 0.0;
	  eval[i] = GauXC::constants::sqrt_pi_ov_2<> * std::erf(sqrt_t) * inv_sqrt_t;
	} else {
	  const double* boys_m = (boys_table + M * DEFAULT_LD_TABLE * DEFAULT_NSEGMENT);
	  constexpr double deltaT = double(DEFAULT_MAX_T) / DEFAULT_NSEGMENT;
	  constexpr double one_over_deltaT = 1 / deltaT;

	  int iseg = std::floor(T[i] * one_over_deltaT);
	  const double* boys_seg = boys_m + iseg * DEFAULT_LD_TABLE;
	  
	  const double ratio = (2 * iseg + 1);
	  const double fact  = 2.0 / deltaT;
	  
	  double xt = fact * T[i] - ratio;
	  double _rec = 1.0;
	  double _val = boys_seg[0];
	  for(int j = 1; j < DEFAULT_NCHEB + 1; ++j) {
	    _rec = _rec * xt;
	    _val += _rec * boys_seg[j];
	  }

	  T_inv_e[i] = 0.5 * std::exp(-T[i]);
	  eval[i] = _val;
	}
      } else {
	const double t_inv = 1./T[i];
	double _val = GauXC::constants::sqrt_pi_ov_2<> * std::sqrt(t_inv);
      
	for(int j = 1; j < M + 1; ++j) {
	  _val *= ((j - 0.5) * t_inv);
	}

	T_inv_e[i] = 0.0;
	eval[i] = _val;
      }
    }
  }
}

// Scalar types
#define SCALAR_TYPE double
#define SCALAR_LENGTH 1

#define SCALAR_SET1(x) (x)

#define SCALAR_LOAD(x) *(x)
#define SCALAR_STORE(x, y) (*(x) = y)

#define SCALAR_ADD(x, y) (x + y)
#define SCALAR_SUB(x, y) (x - y)

#define SCALAR_MUL(x, y) (x * y)
#define SCALAR_FMA(x, y, z) (z + x * y)
#define SCALAR_FNMA(x, y, z) (z - x * y)

#define SCALAR_RECIPROCAL(x) (1.0 / (1.0 * x))

#define SCALAR_DUPLICATE(x) (*(x))

// AVX-512 SIMD Types
#if __AVX512F__ 

  #include <zmmintrin.h>
  
  #define SIMD_TYPE __m512d
  
  #define SIMD_LENGTH 8
  
  #define SIMD_ZERO() _mm512_setzero_pd()
  #define SIMD_SET1(x) _mm512_set1_pd(x)
  
  #define SIMD_ALIGNED_LOAD(x) _mm512_load_pd(x)
  #define SIMD_UNALIGNED_LOAD(x) _mm512_loadu_pd(x)
  
  #define SIMD_ALIGNED_STORE(x, y) _mm512_store_pd(x, y)
  #define SIMD_UNALIGNED_STORE(x, y) _mm512_storeu_pd(x, y)
  
  #define SIMD_ADD(x, y) _mm512_add_pd(x, y)
  #define SIMD_SUB(x, y) _mm512_sub_pd(x, y)
  
  #define SIMD_MUL(x, y) _mm512_mul_pd(x, y)
  #define SIMD_FMA(x, y, z) _mm512_fmadd_pd(x, y, z)
  #define SIMD_FNMA(x, y, z) _mm512_fnmadd_pd(x, y, z)
  
  #define SIMD_DUPLICATE(x) _mm512_broadcast_f64x4(_mm256_broadcast_sd(x))

// AVX-256 SIMD Types
#elif __AVX__ || __AVX2__

  #include <immintrin.h>
  
  #define SIMD_TYPE __m256d
  
  #define SIMD_LENGTH 4
  
  #define SIMD_ZERO() _mm256_setzero_pd()
  #define SIMD_SET1(x) _mm256_set1_pd(x)
  
  #define SIMD_ALIGNED_LOAD(x) _mm256_load_pd(x)
  #define SIMD_UNALIGNED_LOAD(x) _mm256_loadu_pd(x)
  
  #define SIMD_ALIGNED_STORE(x, y) _mm256_store_pd(x, y)
  #define SIMD_UNALIGNED_STORE(x, y) _mm256_storeu_pd(x, y)
  
  #define SIMD_ADD(x, y) _mm256_add_pd(x, y)
  #define SIMD_SUB(x, y) _mm256_sub_pd(x, y)
  
  #define SIMD_MUL(x, y) _mm256_mul_pd(x, y)
  #define SIMD_FMA(x, y, z) _mm256_fmadd_pd(x, y, z)
  #define SIMD_FNMA(x, y, z) _mm256_fnmadd_pd(x, y, z)
  
  #define SIMD_DUPLICATE(x) _mm256_broadcast_sd(x)

// SSE SIMD Types
#elif __SSE__ || __SSE2__ || __SSE3__

  #include <immintrin.h>
  
  #define SIMD_TYPE __m128d
  
  #define SIMD_LENGTH 2
  
  #define SIMD_ZERO() _mm_setzero_pd()
  #define SIMD_SET1(x) _mm_set1_pd(x)
  
  #define SIMD_ALIGNED_LOAD(x) _mm_load_pd(x)
  #define SIMD_UNALIGNED_LOAD(x) _mm_loadu_pd(x)
  
  #define SIMD_ALIGNED_STORE(x, y) _mm_store_pd(x, y)
  #define SIMD_UNALIGNED_STORE(x, y) _mm_storeu_pd(x, y)
  
  #define SIMD_ADD(x, y) _mm_add_pd(x, y)
  #define SIMD_SUB(x, y) _mm_sub_pd(x, y)
  
  #define SIMD_MUL(x, y) _mm_mul_pd(x, y)
  #define SIMD_FMA(x, y, z) _mm_fmadd_pd(x, y, z)
  #define SIMD_FNMA(x, y, z) _mm_fnmadd_pd(x, y, z)
  
  #define SIMD_DUPLICATE(x) _mm_loaddup_pd(x)

// Scalar SIMD Emulation
#else

  #warning "Warning: ISA Not Specified: Using Scalar Code"

  #define SIMD_TYPE double
  
  #define SIMD_LENGTH 1
  
  #define SIMD_ZERO() 0.0
  #define SIMD_SET1(x) x
  
  #define SIMD_ALIGNED_LOAD(x) *(x)
  #define SIMD_UNALIGNED_LOAD(x) *(x)
  
  #define SIMD_ALIGNED_STORE(x, y) *(x) = y
  #define SIMD_UNALIGNED_STORE(x, y) *(x) = y
  
  #define SIMD_ADD(x, y) x + y
  #define SIMD_SUB(x, y) x - y
  
  #define SIMD_MUL(x, y) x * y
  #define SIMD_FMA(x, y, z) z + x * y 
  #define SIMD_FNMA(x, y, z) z - x * y
  
  #define SIMD_DUPLICATE(x) *(x)

#endif

#if 0
#ifdef X86_SCALAR
#elif defined(X86_SSE)
#elif defined(X86_AVX)
#elif defined(X86_AVX512)
#else
  #error "That ISA is not recognized!!!\n"
#endif
#endif

