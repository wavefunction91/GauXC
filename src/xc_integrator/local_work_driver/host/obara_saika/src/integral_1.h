#ifndef __MY_INTEGRAL_1
#define __MY_INTEGRAL_1

#include "integral_1.h"

void integral_1(size_t npts,
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
