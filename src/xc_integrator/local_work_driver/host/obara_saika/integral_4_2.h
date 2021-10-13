#ifdef __MY_INTEGRAL_4_2
#define __MY_INTEGRAL_4_2

void integral_4_2(int npts,
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
