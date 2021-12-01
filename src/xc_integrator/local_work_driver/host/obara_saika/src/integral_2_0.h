#ifndef __MY_INTEGRAL_2_0
#define __MY_INTEGRAL_2_0

#include "integral_2_0.h"

void integral_2_0(size_t npts,
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
