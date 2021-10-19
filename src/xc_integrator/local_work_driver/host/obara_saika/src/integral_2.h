#ifndef __MY_INTEGRAL_2
#define __MY_INTEGRAL_2

#include "integral_2.h"

void integral_2(size_t npts,
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
