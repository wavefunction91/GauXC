#ifdef __MY_INTEGRAL_3_0
#define __MY_INTEGRAL_3_0

void integral_3_0(int npts,
                  shell shellA,
                  shell shellB,
                  point *points,
                  double *Xi,
                  double *Xj,
                  int ldX,
                  double *Gi,
                  double *Gj,
                  int ldG, 
                  double *weights);

#endif
