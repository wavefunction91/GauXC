#ifndef __MY_INTEGRAL_4_0
#define __MY_INTEGRAL_4_0

#include "integral_4_0.h"

void integral_4_0(int npts,
                  shells shellA,
                  shells shellB,
                  point *points,
                  double *Xi,
                  double *Xj,
                  int ldX,
                  double *Gi,
                  double *Gj,
                  int ldG, 
                  double *weights);

#endif
