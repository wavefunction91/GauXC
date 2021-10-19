#include "chebyshev_boys_function.hpp"
#include "boys_computation.h"
#include <cmath>
#include <algorithm>
#include <numeric>

namespace GauXC {

std::unique_ptr<ChebyshevBoysEvaluator> chebyshev_boys_instance;
void gauxc_boys_init(int ncheb, int maxM, int nseg, double minT, double maxT) {
  chebyshev_boys_instance = std::make_unique<ChebyshevBoysEvaluator>(
    ncheb, maxM, nseg, minT, maxT
  );
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



ChebyshevBoysEvaluator::
  ChebyshevBoysEvaluator( int ncheb, int maxM, int nseg, double minT, double maxT ):
  ncheb_(ncheb), maxM_(maxM), nseg_(nseg), min_t_thresh_(minT),
  max_t_thresh_(maxT) { 
  
  table_.resize( (ncheb_+1) * nseg_ * (maxM_+1) );
  ldtable_ = ncheb_+1;

  generate_boys_table( ncheb_, maxM_, max_t_thresh_, nseg_, table_.data(),
    ldtable_ );

}






void monomial_expand( size_t npts, int npoly, const double* coeff, const double *x, 
		      double a, double b, double* eval ) {

  const int n = npoly+1;
  const double sum = a+b;
  const double diff = b-a;
  const double ratio = sum / diff;
  const double fact = 2. / diff;
  for( size_t j = 0; j < npts; ++j ) {
    double xt = fact * x[j] - ratio;
    double xt_use = xt;
    double _val = coeff[0];
    for( int i = 1; i < n; ++i ) {
      _val += xt_use * coeff[i];
      xt_use *= xt;
    }
    eval[j] = _val;
  }

}

void boys_chebyshev( int npts, int m, const double* T, int ncheb, int nseg, double maxT, const double* boys_table, int ld, double* eval ) {
  const double* boys_m = boys_table + m * ld * nseg;

  const double deltaT = maxT / nseg;
  for( int i = 0; i < npts; ++i ) {
    const double tval = T[i];
    if( tval > maxT ) eval[i] = boys_asymp(m,tval);
    else {
      int iseg = std::floor( tval / deltaT);
      const double* boys_seg = boys_m + iseg * ld;

      const double a = iseg * deltaT;
      const double b = a + deltaT;
      monomial_expand( 1, ncheb, boys_seg, T+i, a, b, eval+i );
    }
  }
}

void ChebyshevBoysEvaluator::eval( size_t npts, int m, const double* T, 
  double* FmT ) {
  boys_chebyshev( npts, m, T, ncheb_, nseg_, max_t_thresh_, table_.data(),
    ldtable_, FmT );
}

}
