#ifndef __MY_INTEGRAL_4_4
#define __MY_INTEGRAL_4_4

#include "integral_4_4.h"

void integral_4_4(int npts,
                  shells shellA,
                  shells shellB,
                  point *points,
                  double *Xi,
                  int ldX,
                  double *Gi,
                  int ldG, 
                  double *weights);

#endif
