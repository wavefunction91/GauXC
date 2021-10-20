#ifndef __MY_INTEGRAL_3_1
#define __MY_INTEGRAL_3_1

#include "integral_3_1.h"

void integral_3_1(size_t npts,
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
