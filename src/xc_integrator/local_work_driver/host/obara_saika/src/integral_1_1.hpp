#ifndef __MY_INTEGRAL_1_1
#define __MY_INTEGRAL_1_1

namespace XCPU {
void integral_1_1(size_t npts,
                  shell_pair shpair,
                  double *points,
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