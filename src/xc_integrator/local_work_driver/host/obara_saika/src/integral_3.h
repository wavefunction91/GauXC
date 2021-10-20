#ifndef __MY_INTEGRAL_3
#define __MY_INTEGRAL_3

#include "integral_3.h"

void integral_3(size_t npts,
               shells shellA,
               point *points,
               double *Xi,
               int stX,
               int ldX,
               double *Gi,
               int stG, 
               int ldG, 
               double *weights);

#endif
