#ifndef __MY_INTEGRAL_4
#define __MY_INTEGRAL_4

#include "../include/integral_data_types.hpp"
namespace XCPU {
void integral_4(size_t npts,
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
