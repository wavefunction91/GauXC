#ifndef __MY_INTEGRAL_4_1
#define __MY_INTEGRAL_4_1

#include "integral_4_1.h"

void integral_4_1(int npts,
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
