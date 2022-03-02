#ifndef __MY_INTEGRAL_4_0
#define __MY_INTEGRAL_4_0

#include "../include/integral_data_types.hpp"
namespace XCPU {
void integral_4_0(size_t npts,
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
