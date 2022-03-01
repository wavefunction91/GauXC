#ifndef __MY_INTEGRAL_1
#define __MY_INTEGRAL_1

namespace XCPU {
void integral_1(size_t npts,
               double *points,
               point rA,
               point rB,
               int nprim_pairs,
               prim_pair *prim_pairs,
               double *Xi,
               int ldX,
               double *Gi,
               int ldG, 
               double *weights, 
               double *boys_table);
}

#endif
