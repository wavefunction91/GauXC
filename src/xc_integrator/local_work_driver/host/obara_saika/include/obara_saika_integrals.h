#ifndef __MY_INTEGRAL_OBARA_SAIKA
#define __MY_INTEGRAL_OBARA_SAIKA

void compute_integral_shell_pair(int npts,
                  shells shellA,
                  shells shellB,
                  point *points,
                  double *Xi,
                  double *Xj,
                  int ldX,
                  double *Gi,
                  double *Gj,
                  int ldG, 
                  double *weights);

#endif
