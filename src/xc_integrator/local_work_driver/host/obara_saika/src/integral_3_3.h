#ifndef __MY_INTEGRAL_3_3
#define __MY_INTEGRAL_3_3

#include "integral_3_3.h"

void integral_3_3(size_t npts,
                  shell_pair shpair,
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