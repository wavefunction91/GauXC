#ifndef __MY_INTEGRAL_4_2
#define __MY_INTEGRAL_4_2

#include "integral_4_2.h"

void integral_4_2(size_t npts,
                  shells shellA,
                  shells shellB,
                  point *points,
                  double *Xi,
                  double *Xj,
                  int stX,
                  int ldX,
                  double *Gi,
                  double *Gj,
                  int stG, 
                  int ldG, 
                  double *weights);

#endif
