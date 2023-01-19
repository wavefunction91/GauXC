#ifndef __MY_INTEGRAL_1_0
#define __MY_INTEGRAL_1_0

#include "../include/cpu/integral_data_types.hpp"
namespace XCPU {
void integral_1_0(size_t npts,
                  double *points,
                  point rA,
                  point rB,
                  int nprim_pairs,
                  prim_pair *prim_pairs,
                  double *Xi,
                  double *Xj,
                  int ldX,
                  double *Gi,
                  double *Gj,
                  int ldG, 
                  double *weights, 
                  double *boys_table);
}

#endif
