#ifndef __MY_INTEGRAL_OBARA_SAIKA
#define __MY_INTEGRAL_OBARA_SAIKA

namespace XCPU {
void generate_shell_pair( const shells& A, const shells& B, prim_pair *prim_pairs);
void compute_integral_shell_pair(int is_diag,
                  size_t npts,
                  double *points,
                  int lA,
                  int lB,
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
