#ifndef __MY_INTEGRAL_OBARA_SAIKA
#define __MY_INTEGRAL_OBARA_SAIKA
namespace XGPU {
void generate_shell_pair( const shells& A, const shells& B, shell_pair& AB);
void compute_integral_shell_pair(size_t npts,
                             int is_diag,
                             int lA,
                             int lB,
                             point rA,
                             point rB,
                             point rAB,
                             int nprim_pair,
                             prim_pair *ppair,
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
