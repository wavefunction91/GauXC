#include "chebyshev_boys_function.hpp"
#include "boys_computation.h"
#include <cmath>
#include <algorithm>
#include <numeric>

#define MIN(a,b)			\
  ({ __typeof__ (a) _a = (a);	        \
  __typeof__ (b) _b = (b);		\
  _a < _b ? _a : _b; })

namespace GauXC {

std::unique_ptr<detail::default_chebyshev_type> chebyshev_boys_instance;
void gauxc_boys_init() {
  chebyshev_boys_instance = std::make_unique<detail::default_chebyshev_type>();
}
void gauxc_boys_finalize() {
  chebyshev_boys_instance.reset();
}












template <typename Op>
void cheby_coeff( int ncheb, const Op& f, double a, double b, double* c ) {

  const int n = ncheb+1;
  const double pi_ov_2n = M_PI / (2 * n);

  // Generate function table
  std::vector<double> f_table(n);
  for( int i = 0; i < n; ++i ) {
    double x = std::cos( (2.*(i+1)-1) * pi_ov_2n ); // Chebyshev node
    x = 0.5 * ( a+b + (b-a)*x );
    f_table[i] = f(x);
  }

  c[0] = std::accumulate( f_table.begin(), f_table.end(),0. ) / n;
  for( int i = 1; i < n; ++i ) {
    double _val = 0.;
    for( int j = 0; j < n; ++j ) {
      // f(x_j) * T_i(x_j)
      _val += f_table[j] * std::cos( i * (2*(j+1)-1) * pi_ov_2n );
    }
    c[i] = 2.0 * _val / n;
  }

}

void cheby_to_monomial_coeff( int ncheb, double *coeff ) {

  const int n = ncheb+1;
  int64_t i_fact = 1; // i!
  int64_t t_fact = 1; // 2^(i-1)
  for(int i = 0; i < n; ++i) {
    if(i)     i_fact *= i; // Update i!
    if(i > 1) t_fact *= 2; // Update 2^(i-1)

    double _val = 0;
    if(!i) {
      int m1_fact = 1; // (-1)^(j/2)
      for( int j = i; j < n; j += 2 ) {
        _val += m1_fact * coeff[j];
	      m1_fact *= -1; // Update (-1)^(j/2)
      }
    } else {

      // OEIS A008310
      // a(i,j) = 2^(i-1) * j * (-1)^((j-i)/2) * ((i+j)/2-1)! / (i! * ((j-i)/2)!)
      // Coeff tform is UT, able to be performed in place
      int m1_term = 1;
      for( int j = i; j < n; j += 2 ) {
	      const int f_up = (i+j)/2 - 1;
        const int f_lo = (j-i)/2;
	      // This comptes the ratio of factorials as a 
	      // Pochhammer falling factorial
	      int f_term = 1;
	      for( int k = f_lo+1; k <= f_up; ++k ) f_term *= k;
        _val += t_fact * j * m1_term * double(f_term) / double(i_fact) * coeff[j];
	      m1_term *= -1;
      }

    }

    coeff[i] = _val;
  }

}

void generate_boys_table(int ncheb, int maxM, double maxT, int nseg, 
  double* cheb_coeff_table, int ld) {
  // cheb_coeff_table is (ld, nseg, maxM+1) ld >= (ncheb+1)

  const double deltaT = maxT / nseg;
  for( int m = 0; m <= maxM; ++m ) {
    double* coeff_m = cheb_coeff_table + m * ld * nseg; // table offset for current m
    for( int iseg = 0; iseg < nseg; ++iseg ) {
      double* coeff_seg = coeff_m + iseg * ld;

      const double a = iseg * deltaT;
      const double b = a + deltaT;

      auto f = [=](double x){ return boys_reference(m,x); };
      cheby_coeff( ncheb, f, a, b, coeff_seg ); // Generate coeff in Chebyshev basis
      cheby_to_monomial_coeff( ncheb, coeff_seg );   // Convert to monomial basis
    }
  }

}


template <uint32_t NCheb, uint32_t MaxM, uint32_t MaxT, uint32_t NSegment, 
          uint32_t LDTable>
ChebyshevBoysEvaluator<NCheb,MaxM,MaxT,NSegment,LDTable>::ChebyshevBoysEvaluator() {
  
  table_.resize( (NCheb+1) * NSegment * (MaxM+1) );

  generate_boys_table( NCheb, MaxM, MaxT, NSegment, table_.data(),
    LDTable );

}





template <uint32_t NPoly>
void monomial_expand( size_t npts, const double* coeff, const double *x, 
		      double a, double b, double* eval ) {

  constexpr int n = NPoly+1;
  const double sum = a+b;
  const double diff = b-a;
  const double ratio = sum / diff;
  const double fact = 2. / diff;

  double xp[n]; xp[0] = 1.;

  for( size_t j = 0; j < npts; ++j ) {
    double xt = fact * x[j] - ratio;

    #pragma unroll
    for( int i = 1; i < n; ++i ) xp[i] = xp[i-1] * xt;

    double _val = 0.;
    #pragma unroll
    for( int i = 0; i < n; ++i ) _val += xp[i] * coeff[i];

    eval[j] = _val;
  }

}

template <uint32_t NCheb, uint32_t MaxM, uint32_t MaxT, uint32_t NSegment, 
          uint32_t LDTable>
void boys_chebyshev( int npts, int m, const double* T, const double* boys_table, double* eval ) {
  const double* boys_m = boys_table + m * LDTable * NSegment;

  constexpr double deltaT = double(MaxT) / NSegment;
#if 1
  for( int i = 0; i < npts; ++i ) {
    const double tval = T[i];
    if( tval > MaxT ) {eval[i] = boys_asymp(m,tval); }
    else {
      int iseg = std::floor( tval / deltaT);
      const double* boys_seg = boys_m + iseg * LDTable;

      const double a = iseg * deltaT;
      const double b = a + deltaT;
      monomial_expand<NCheb>( 1, boys_seg, T+i, a, b, eval+i );
    }
  }
#else
  //double temp[NPTS_LOCAL];
  for( int i_st = 0; i_st < npts; i_st += NPTS_LOCAL ) {
    int ndo = MIN( NPTS_LOCAL, npts - i_st );
    // Get min/max element
    auto [min_t_it,max_t_it] = std::minmax_element( T + i_st, T + i_st + ndo );

    if( *min_t_it >= MaxT ) {
      // All asymptotic case
      //boys_asymp( ndo, m, T + i_st, eval + i_st );
    } else {
      // Mixed case - need to figure out a way to vectorize Monomial expansion
      for( int i = 0; i < npts; ++i ) {
        const double tval = T[i];
        if( tval > MaxT ) {eval[i] = boys_asymp(m,tval);}
        else {
          int iseg = std::floor( tval / deltaT);
          const double* boys_seg = boys_m + iseg * LDTable;
          //printf("%d\n",iseg);

          const double a = iseg * deltaT;
          const double b = a + deltaT;
          monomial_expand<NCheb>( 1, boys_seg, T+i, a, b, eval+i );
        }
      }
    }
  }
#endif
}

template <uint32_t NCheb, uint32_t MaxM, uint32_t MaxT, uint32_t NSegment, 
          uint32_t LDTable>
void ChebyshevBoysEvaluator<NCheb,MaxM,MaxT,NSegment,LDTable>::eval( size_t npts, int m, const double* T, 
  double* FmT ) {
  boys_chebyshev<NCheb,MaxM,MaxT,NSegment,LDTable>( npts, m, T, table_.data(), FmT );
}

template class ChebyshevBoysEvaluator< detail::default_ncheb, detail::default_max_m, detail::default_max_t >;

}
