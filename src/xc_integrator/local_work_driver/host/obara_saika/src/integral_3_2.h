#ifndef __MY_INTEGRAL_3_2
#define __MY_INTEGRAL_3_2

#include "integral_3_2.h"

void integral_3_2(size_t npts,
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
