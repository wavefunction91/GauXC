#ifndef __MY_INTEGRAL_3_3
#define __MY_INTEGRAL_3_3

#include "integral_3_3.h"

void integral_3_3(int npts,
                  shells shellA,
                  shells shellB,
                  point *points,
                  double *Xi,
                  int ldX,
                  double *Gi,
                  int ldG, 
                  double *weights);

#endif
