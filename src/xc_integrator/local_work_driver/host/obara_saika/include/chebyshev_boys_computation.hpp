#ifndef _MY_CHEBYSHEV_BOYS_FUNCTION
#define _MY_CHEBYSHEV_BOYS_FUNCTION

#include <iostream>

#define DEFAULT_NCHEB  7
#define DEFAULT_MAX_M 16
#define DEFAULT_MAX_T 30

#define DEFAULT_NSEGMENT ((DEFAULT_MAX_T * DEFAULT_NCHEB) / 2)
#define DEFAULT_LD_TABLE (DEFAULT_NCHEB + 1)

namespace GauXC {

  // create tables
  
  void gauxc_boys_init();
  
  void gauxc_boys_finalize();

  // get values

  template <int M>
  double gauxc_boys_element(double T);

  template <int M>
  void gauxc_boys_elements(size_t npts, double *T, double *eval);

#define BOYS_FUNCTION_HEADER(M)						\
  extern template double gauxc_boys_element<M>(double);			\
  extern template void gauxc_boys_elements<M>(size_t, double*, double*);

  BOYS_FUNCTION_HEADER( 0);
  BOYS_FUNCTION_HEADER( 1);
  BOYS_FUNCTION_HEADER( 2);
  BOYS_FUNCTION_HEADER( 3);
  BOYS_FUNCTION_HEADER( 4);
  BOYS_FUNCTION_HEADER( 5);
  BOYS_FUNCTION_HEADER( 6);
  BOYS_FUNCTION_HEADER( 7);
  BOYS_FUNCTION_HEADER( 8);
  BOYS_FUNCTION_HEADER( 9);
  BOYS_FUNCTION_HEADER(10);
  BOYS_FUNCTION_HEADER(11);
  BOYS_FUNCTION_HEADER(12);
  BOYS_FUNCTION_HEADER(13);
  BOYS_FUNCTION_HEADER(14);
  BOYS_FUNCTION_HEADER(15);
}

#endif
