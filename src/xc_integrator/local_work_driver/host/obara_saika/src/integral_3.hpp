#ifndef __MY_INTEGRAL_3
#define __MY_INTEGRAL_3

#include "../include/cpu/integral_data_types.hpp"
namespace XCPU {
void integral_3(size_t npts,
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
