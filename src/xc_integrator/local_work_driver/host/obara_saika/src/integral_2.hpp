#ifndef __MY_INTEGRAL_2
#define __MY_INTEGRAL_2

namespace XCPU {
void integral_2(size_t npts,
               shell_pair shpair,
               double *points,
               double *Xi,
               int ldX,
               double *Gi,
               int ldG, 
               double *weights, 
               double *boys_table);
}

#endif