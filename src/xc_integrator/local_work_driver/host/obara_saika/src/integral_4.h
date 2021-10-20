#ifndef __MY_INTEGRAL_4
#define __MY_INTEGRAL_4

#include "integral_4.h"

void integral_4(size_t npts,
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
