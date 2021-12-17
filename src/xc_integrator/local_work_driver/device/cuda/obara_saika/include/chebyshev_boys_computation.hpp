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
  
  double *gauxc_boys_init();
  
  void gauxc_boys_finalize();

}

#endif
