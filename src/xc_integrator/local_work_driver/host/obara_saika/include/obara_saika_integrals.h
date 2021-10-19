#ifndef __MY_INTEGRAL_OBARA_SAIKA
#define __MY_INTEGRAL_OBARA_SAIKA

void compute_integral_shell_pair(size_t npts,
                  int i,
                  int j,
                  shells *shell_list,
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
