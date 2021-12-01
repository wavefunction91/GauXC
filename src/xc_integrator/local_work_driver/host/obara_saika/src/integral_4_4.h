#ifndef __MY_INTEGRAL_4_4
#define __MY_INTEGRAL_4_4

#include "integral_4_4.h"

void integral_4_4(size_t npts,
                  shell_pair shpair,
                  double *points,
                  double *Xi,
                  double *Xj,
                  int ldX,
                  double *Gi,
                  double *Gj,
                  int ldG, 
                  double *weights);

#endif
