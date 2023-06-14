#ifndef __RYS_INTEGRALS
#define __RYS_INTEGRALS

typedef struct {
  double x, y, z;
} point;

typedef struct {
  double alpha, coeff;
} coefficients;

typedef struct {
  point origin;
  coefficients *coeff;
  int m, L;
} shells;

typedef struct {
  point P;
  point PA;
  point PB;

  double K;
  double gamma;
  double coeff_prod;
} prim_pair;

typedef struct {
  int lA;
  int lB;
  int nprim_pair;
  point rAB;
  prim_pair* prim_pairs;
} shell_pair;

#ifdef __cplusplus
extern "C" {
#endif
void compute_integral(int n, shells *shell_list, int m, point *points, double *output);
void compute_integral_shell_pair( int npts, shells sh0, shells sh1, 
                                  point *points, double* matrix ); 
void compute_integral_shell_pair_pre( int npts, shell_pair shpair,
                                      point* points, double* matrix );
#ifdef __cplusplus
}
#endif

#endif
