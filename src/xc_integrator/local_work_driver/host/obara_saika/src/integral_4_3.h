#ifndef __MY_INTEGRAL_4_3
#define __MY_INTEGRAL_4_3

#include "integral_4_3.h"

void integral_4_3(size_t npts,
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
