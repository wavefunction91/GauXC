#include <math.h>
#include "boys_computation.h"
#include "integral_data_types.h"
#include "config_obara_saika.h"

#define PI 3.14159265358979323846

#define MIN(a,b)			\
  ({ __typeof__ (a) _a = (a);	        \
  __typeof__ (b) _b = (b);		\
  _a < _b ? _a : _b; })

void integral_2(size_t npts,
               shell_pair shpair,
               double *_points,
               double *Xi,
               int ldX,
               double *Gi,
               int ldG, 
               double *weights) {
   __attribute__((__aligned__(64))) double buffer[31 * NPTS_LOCAL + 5 * NPTS_LOCAL + NPTS_LOCAL];

   double *temp = (buffer + 0);
   double *FmT = (buffer + 31 * NPTS_LOCAL);
   double *Tval = (buffer + 31 * NPTS_LOCAL + 5 * NPTS_LOCAL);

   size_t npts_upper = NPTS_LOCAL * (npts / NPTS_LOCAL);
   size_t p_outer = 0;
   for(p_outer = 0; p_outer < npts_upper; p_outer += NPTS_LOCAL) {
      double *_point_outer = (_points + p_outer);

      double xA = shpair.rA.x;
      double yA = shpair.rA.y;
      double zA = shpair.rA.z;

      for(int i = 0; i < 31 * NPTS_LOCAL; i += SIMD_LENGTH) SIMD_ALIGNED_STORE((temp + i), SIMD_ZERO());

      for(int ij = 0; ij < shpair.nprim_pair; ++ij) {
         double RHO = shpair.prim_pairs[ij].gamma;
         double RHO_INV = 1.0 / RHO;

         constexpr double X_PA = 0.0;
         constexpr double Y_PA = 0.0;
         constexpr double Z_PA = 0.0;

         double eval = shpair.prim_pairs[ij].coeff_prod * shpair.prim_pairs[ij].K;

         // Evaluate T Values
         for(size_t p_inner = 0; p_inner < NPTS_LOCAL; p_inner += SIMD_LENGTH) {
            SIMD_TYPE xC = SIMD_UNALIGNED_LOAD((_point_outer + p_inner + 0 * npts));
            SIMD_TYPE yC = SIMD_UNALIGNED_LOAD((_point_outer + p_inner + 1 * npts));
            SIMD_TYPE zC = SIMD_UNALIGNED_LOAD((_point_outer + p_inner + 2 * npts));

            SIMD_TYPE X_PC = SIMD_SUB(SIMD_DUPLICATE(&(xA)), xC);
            SIMD_TYPE Y_PC = SIMD_SUB(SIMD_DUPLICATE(&(yA)), yC);
            SIMD_TYPE Z_PC = SIMD_SUB(SIMD_DUPLICATE(&(zA)), zC);

            X_PC = SIMD_MUL(X_PC, X_PC);
            X_PC = SIMD_FMA(Y_PC, Y_PC, X_PC);
            X_PC = SIMD_FMA(Z_PC, Z_PC, X_PC);
            X_PC = SIMD_MUL(SIMD_DUPLICATE(&(RHO)), X_PC);
            SIMD_ALIGNED_STORE((Tval + p_inner), X_PC);
         }

         // Evaluate Boys function
         GauXC::boys_function<0>(NPTS_LOCAL, Tval, FmT + 0 * NPTS_LOCAL);
         GauXC::boys_function<1>(NPTS_LOCAL, Tval, FmT + 1 * NPTS_LOCAL);
         GauXC::boys_function<2>(NPTS_LOCAL, Tval, FmT + 2 * NPTS_LOCAL);
         GauXC::boys_function<3>(NPTS_LOCAL, Tval, FmT + 3 * NPTS_LOCAL);
         GauXC::boys_function<4>(NPTS_LOCAL, Tval, FmT + 4 * NPTS_LOCAL);

         // Evaluate VRR Buffer
         for(size_t p_inner = 0; p_inner < NPTS_LOCAL; p_inner += SIMD_LENGTH) {
            SIMD_TYPE xC = SIMD_UNALIGNED_LOAD((_point_outer + p_inner + 0 * npts));
            SIMD_TYPE yC = SIMD_UNALIGNED_LOAD((_point_outer + p_inner + 1 * npts));
            SIMD_TYPE zC = SIMD_UNALIGNED_LOAD((_point_outer + p_inner + 2 * npts));

            SIMD_TYPE X_PC = SIMD_SUB(SIMD_DUPLICATE(&(xA)), xC);
            SIMD_TYPE Y_PC = SIMD_SUB(SIMD_DUPLICATE(&(yA)), yC);
            SIMD_TYPE Z_PC = SIMD_SUB(SIMD_DUPLICATE(&(zA)), zC);

            SIMD_TYPE tx, ty, t00, t01, t02, t03, t04, t10, t11, t12, t13, t20, t21, t22, t30, t31, t40;

            t00 = SIMD_ALIGNED_LOAD((FmT + p_inner + 0 * NPTS_LOCAL));
            t00 = SIMD_MUL(SIMD_DUPLICATE(&(eval)), t00);
            t01 = SIMD_ALIGNED_LOAD((FmT + p_inner + 1 * NPTS_LOCAL));
            t01 = SIMD_MUL(SIMD_DUPLICATE(&(eval)), t01);
            t02 = SIMD_ALIGNED_LOAD((FmT + p_inner + 2 * NPTS_LOCAL));
            t02 = SIMD_MUL(SIMD_DUPLICATE(&(eval)), t02);
            t03 = SIMD_ALIGNED_LOAD((FmT + p_inner + 3 * NPTS_LOCAL));
            t03 = SIMD_MUL(SIMD_DUPLICATE(&(eval)), t03);
            t04 = SIMD_ALIGNED_LOAD((FmT + p_inner + 4 * NPTS_LOCAL));
            t04 = SIMD_MUL(SIMD_DUPLICATE(&(eval)), t04);
            t10 = SIMD_MUL(SIMD_DUPLICATE(&(X_PA)), t00);
            t10 = SIMD_FNMA(X_PC, t01, t10);
            t11 = SIMD_MUL(SIMD_DUPLICATE(&(X_PA)), t01);
            t11 = SIMD_FNMA(X_PC, t02, t11);
            t12 = SIMD_MUL(SIMD_DUPLICATE(&(X_PA)), t02);
            t12 = SIMD_FNMA(X_PC, t03, t12);
            t13 = SIMD_MUL(SIMD_DUPLICATE(&(X_PA)), t03);
            t13 = SIMD_FNMA(X_PC, t04, t13);
            t20 = SIMD_MUL(SIMD_DUPLICATE(&(X_PA)), t10);
            t20 = SIMD_FNMA(X_PC, t11, t20);
            tx = SIMD_SUB(t00, t01);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t20 = SIMD_FMA(tx, ty, t20);
            t21 = SIMD_MUL(SIMD_DUPLICATE(&(X_PA)), t11);
            t21 = SIMD_FNMA(X_PC, t12, t21);
            tx = SIMD_SUB(t01, t02);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t21 = SIMD_FMA(tx, ty, t21);
            t22 = SIMD_MUL(SIMD_DUPLICATE(&(X_PA)), t12);
            t22 = SIMD_FNMA(X_PC, t13, t22);
            tx = SIMD_SUB(t02, t03);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t22 = SIMD_FMA(tx, ty, t22);
            tx = SIMD_ALIGNED_LOAD((temp + 0 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t20);
            SIMD_ALIGNED_STORE((temp + 0 * NPTS_LOCAL + p_inner), tx);
            t30 = SIMD_MUL(SIMD_DUPLICATE(&(X_PA)), t20);
            t30 = SIMD_FNMA(X_PC, t21, t30);
            tx = SIMD_SUB(t10, t11);
            ty = SIMD_SET1(0.5 * 2);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t30 = SIMD_FMA(tx, ty, t30);
            t31 = SIMD_MUL(SIMD_DUPLICATE(&(X_PA)), t21);
            t31 = SIMD_FNMA(X_PC, t22, t31);
            tx = SIMD_SUB(t11, t12);
            ty = SIMD_SET1(0.5 * 2);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t31 = SIMD_FMA(tx, ty, t31);
            tx = SIMD_ALIGNED_LOAD((temp + 6 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t30);
            SIMD_ALIGNED_STORE((temp + 6 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(X_PA)), t30);
            t40 = SIMD_FNMA(X_PC, t31, t40);
            tx = SIMD_SUB(t20, t21);
            ty = SIMD_SET1(0.5 * 3);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t40 = SIMD_FMA(tx, ty, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 16 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 16 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t30);
            t40 = SIMD_FNMA(Y_PC, t31, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 17 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 17 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t30);
            t40 = SIMD_FNMA(Z_PC, t31, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 18 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 18 * NPTS_LOCAL + p_inner), tx);
            t30 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t20);
            t30 = SIMD_FNMA(Y_PC, t21, t30);
            t31 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t21);
            t31 = SIMD_FNMA(Y_PC, t22, t31);
            tx = SIMD_ALIGNED_LOAD((temp + 7 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t30);
            SIMD_ALIGNED_STORE((temp + 7 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t30);
            t40 = SIMD_FNMA(Y_PC, t31, t40);
            tx = SIMD_SUB(t20, t21);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t40 = SIMD_FMA(tx, ty, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 19 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 19 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t30);
            t40 = SIMD_FNMA(Z_PC, t31, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 20 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 20 * NPTS_LOCAL + p_inner), tx);
            t30 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t20);
            t30 = SIMD_FNMA(Z_PC, t21, t30);
            t31 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t21);
            t31 = SIMD_FNMA(Z_PC, t22, t31);
            tx = SIMD_ALIGNED_LOAD((temp + 8 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t30);
            SIMD_ALIGNED_STORE((temp + 8 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t30);
            t40 = SIMD_FNMA(Z_PC, t31, t40);
            tx = SIMD_SUB(t20, t21);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t40 = SIMD_FMA(tx, ty, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 21 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 21 * NPTS_LOCAL + p_inner), tx);
            t20 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t10);
            t20 = SIMD_FNMA(Y_PC, t11, t20);
            t21 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t11);
            t21 = SIMD_FNMA(Y_PC, t12, t21);
            t22 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t12);
            t22 = SIMD_FNMA(Y_PC, t13, t22);
            tx = SIMD_ALIGNED_LOAD((temp + 1 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t20);
            SIMD_ALIGNED_STORE((temp + 1 * NPTS_LOCAL + p_inner), tx);
            t30 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t20);
            t30 = SIMD_FNMA(Y_PC, t21, t30);
            tx = SIMD_SUB(t10, t11);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t30 = SIMD_FMA(tx, ty, t30);
            t31 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t21);
            t31 = SIMD_FNMA(Y_PC, t22, t31);
            tx = SIMD_SUB(t11, t12);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t31 = SIMD_FMA(tx, ty, t31);
            tx = SIMD_ALIGNED_LOAD((temp + 9 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t30);
            SIMD_ALIGNED_STORE((temp + 9 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t30);
            t40 = SIMD_FNMA(Y_PC, t31, t40);
            tx = SIMD_SUB(t20, t21);
            ty = SIMD_SET1(0.5 * 2);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t40 = SIMD_FMA(tx, ty, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 22 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 22 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t30);
            t40 = SIMD_FNMA(Z_PC, t31, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 23 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 23 * NPTS_LOCAL + p_inner), tx);
            t30 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t20);
            t30 = SIMD_FNMA(Z_PC, t21, t30);
            t31 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t21);
            t31 = SIMD_FNMA(Z_PC, t22, t31);
            tx = SIMD_ALIGNED_LOAD((temp + 10 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t30);
            SIMD_ALIGNED_STORE((temp + 10 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t30);
            t40 = SIMD_FNMA(Z_PC, t31, t40);
            tx = SIMD_SUB(t20, t21);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t40 = SIMD_FMA(tx, ty, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 24 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 24 * NPTS_LOCAL + p_inner), tx);
            t20 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t10);
            t20 = SIMD_FNMA(Z_PC, t11, t20);
            t21 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t11);
            t21 = SIMD_FNMA(Z_PC, t12, t21);
            t22 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t12);
            t22 = SIMD_FNMA(Z_PC, t13, t22);
            tx = SIMD_ALIGNED_LOAD((temp + 2 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t20);
            SIMD_ALIGNED_STORE((temp + 2 * NPTS_LOCAL + p_inner), tx);
            t30 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t20);
            t30 = SIMD_FNMA(Z_PC, t21, t30);
            tx = SIMD_SUB(t10, t11);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t30 = SIMD_FMA(tx, ty, t30);
            t31 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t21);
            t31 = SIMD_FNMA(Z_PC, t22, t31);
            tx = SIMD_SUB(t11, t12);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t31 = SIMD_FMA(tx, ty, t31);
            tx = SIMD_ALIGNED_LOAD((temp + 11 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t30);
            SIMD_ALIGNED_STORE((temp + 11 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t30);
            t40 = SIMD_FNMA(Z_PC, t31, t40);
            tx = SIMD_SUB(t20, t21);
            ty = SIMD_SET1(0.5 * 2);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t40 = SIMD_FMA(tx, ty, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 25 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 25 * NPTS_LOCAL + p_inner), tx);
            t10 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t00);
            t10 = SIMD_FNMA(Y_PC, t01, t10);
            t11 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t01);
            t11 = SIMD_FNMA(Y_PC, t02, t11);
            t12 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t02);
            t12 = SIMD_FNMA(Y_PC, t03, t12);
            t13 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t03);
            t13 = SIMD_FNMA(Y_PC, t04, t13);
            t20 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t10);
            t20 = SIMD_FNMA(Y_PC, t11, t20);
            tx = SIMD_SUB(t00, t01);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t20 = SIMD_FMA(tx, ty, t20);
            t21 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t11);
            t21 = SIMD_FNMA(Y_PC, t12, t21);
            tx = SIMD_SUB(t01, t02);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t21 = SIMD_FMA(tx, ty, t21);
            t22 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t12);
            t22 = SIMD_FNMA(Y_PC, t13, t22);
            tx = SIMD_SUB(t02, t03);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t22 = SIMD_FMA(tx, ty, t22);
            tx = SIMD_ALIGNED_LOAD((temp + 3 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t20);
            SIMD_ALIGNED_STORE((temp + 3 * NPTS_LOCAL + p_inner), tx);
            t30 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t20);
            t30 = SIMD_FNMA(Y_PC, t21, t30);
            tx = SIMD_SUB(t10, t11);
            ty = SIMD_SET1(0.5 * 2);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t30 = SIMD_FMA(tx, ty, t30);
            t31 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t21);
            t31 = SIMD_FNMA(Y_PC, t22, t31);
            tx = SIMD_SUB(t11, t12);
            ty = SIMD_SET1(0.5 * 2);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t31 = SIMD_FMA(tx, ty, t31);
            tx = SIMD_ALIGNED_LOAD((temp + 12 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t30);
            SIMD_ALIGNED_STORE((temp + 12 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t30);
            t40 = SIMD_FNMA(Y_PC, t31, t40);
            tx = SIMD_SUB(t20, t21);
            ty = SIMD_SET1(0.5 * 3);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t40 = SIMD_FMA(tx, ty, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 26 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 26 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t30);
            t40 = SIMD_FNMA(Z_PC, t31, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 27 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 27 * NPTS_LOCAL + p_inner), tx);
            t30 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t20);
            t30 = SIMD_FNMA(Z_PC, t21, t30);
            t31 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t21);
            t31 = SIMD_FNMA(Z_PC, t22, t31);
            tx = SIMD_ALIGNED_LOAD((temp + 13 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t30);
            SIMD_ALIGNED_STORE((temp + 13 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t30);
            t40 = SIMD_FNMA(Z_PC, t31, t40);
            tx = SIMD_SUB(t20, t21);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t40 = SIMD_FMA(tx, ty, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 28 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 28 * NPTS_LOCAL + p_inner), tx);
            t20 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t10);
            t20 = SIMD_FNMA(Z_PC, t11, t20);
            t21 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t11);
            t21 = SIMD_FNMA(Z_PC, t12, t21);
            t22 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t12);
            t22 = SIMD_FNMA(Z_PC, t13, t22);
            tx = SIMD_ALIGNED_LOAD((temp + 4 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t20);
            SIMD_ALIGNED_STORE((temp + 4 * NPTS_LOCAL + p_inner), tx);
            t30 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t20);
            t30 = SIMD_FNMA(Z_PC, t21, t30);
            tx = SIMD_SUB(t10, t11);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t30 = SIMD_FMA(tx, ty, t30);
            t31 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t21);
            t31 = SIMD_FNMA(Z_PC, t22, t31);
            tx = SIMD_SUB(t11, t12);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t31 = SIMD_FMA(tx, ty, t31);
            tx = SIMD_ALIGNED_LOAD((temp + 14 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t30);
            SIMD_ALIGNED_STORE((temp + 14 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t30);
            t40 = SIMD_FNMA(Z_PC, t31, t40);
            tx = SIMD_SUB(t20, t21);
            ty = SIMD_SET1(0.5 * 2);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t40 = SIMD_FMA(tx, ty, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 29 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 29 * NPTS_LOCAL + p_inner), tx);
            t10 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t00);
            t10 = SIMD_FNMA(Z_PC, t01, t10);
            t11 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t01);
            t11 = SIMD_FNMA(Z_PC, t02, t11);
            t12 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t02);
            t12 = SIMD_FNMA(Z_PC, t03, t12);
            t13 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t03);
            t13 = SIMD_FNMA(Z_PC, t04, t13);
            t20 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t10);
            t20 = SIMD_FNMA(Z_PC, t11, t20);
            tx = SIMD_SUB(t00, t01);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t20 = SIMD_FMA(tx, ty, t20);
            t21 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t11);
            t21 = SIMD_FNMA(Z_PC, t12, t21);
            tx = SIMD_SUB(t01, t02);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t21 = SIMD_FMA(tx, ty, t21);
            t22 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t12);
            t22 = SIMD_FNMA(Z_PC, t13, t22);
            tx = SIMD_SUB(t02, t03);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t22 = SIMD_FMA(tx, ty, t22);
            tx = SIMD_ALIGNED_LOAD((temp + 5 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t20);
            SIMD_ALIGNED_STORE((temp + 5 * NPTS_LOCAL + p_inner), tx);
            t30 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t20);
            t30 = SIMD_FNMA(Z_PC, t21, t30);
            tx = SIMD_SUB(t10, t11);
            ty = SIMD_SET1(0.5 * 2);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t30 = SIMD_FMA(tx, ty, t30);
            t31 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t21);
            t31 = SIMD_FNMA(Z_PC, t22, t31);
            tx = SIMD_SUB(t11, t12);
            ty = SIMD_SET1(0.5 * 2);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t31 = SIMD_FMA(tx, ty, t31);
            tx = SIMD_ALIGNED_LOAD((temp + 15 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t30);
            SIMD_ALIGNED_STORE((temp + 15 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t30);
            t40 = SIMD_FNMA(Z_PC, t31, t40);
            tx = SIMD_SUB(t20, t21);
            ty = SIMD_SET1(0.5 * 3);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t40 = SIMD_FMA(tx, ty, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 30 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 30 * NPTS_LOCAL + p_inner), tx);
         }
      }

      for(size_t p_inner = 0; p_inner < NPTS_LOCAL; p_inner += SIMD_LENGTH) {
         double *Xik = (Xi + p_outer + p_inner);
         double *Gik = (Gi + p_outer + p_inner);

         SIMD_TYPE tx, wg, xik, gik;
         tx  = SIMD_ALIGNED_LOAD((temp + 16 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 0 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 0 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 0 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 17 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 0 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 1 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 1 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 18 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 0 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 2 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 2 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 19 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 0 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 3 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 3 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 20 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 0 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 4 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 4 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 21 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 0 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 5 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 5 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 17 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 1 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 0 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 0 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 19 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 1 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 1 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 1 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 20 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 1 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 2 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 2 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 22 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 1 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 3 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 3 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 23 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 1 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 4 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 4 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 24 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 1 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 5 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 5 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 18 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 2 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 0 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 0 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 20 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 2 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 1 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 1 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 21 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 2 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 2 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 2 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 23 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 2 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 3 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 3 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 24 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 2 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 4 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 4 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 25 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 2 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 5 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 5 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 19 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 3 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 0 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 0 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 22 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 3 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 1 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 1 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 23 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 3 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 2 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 2 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 26 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 3 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 3 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 3 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 27 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 3 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 4 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 4 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 28 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 3 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 5 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 5 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 20 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 4 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 0 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 0 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 23 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 4 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 1 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 1 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 24 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 4 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 2 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 2 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 27 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 4 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 3 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 3 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 28 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 4 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 4 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 4 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 29 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 4 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 5 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 5 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 21 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 5 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 0 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 0 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 24 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 5 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 1 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 1 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 25 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 5 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 2 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 2 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 28 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 5 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 3 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 3 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 29 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 5 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 4 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 4 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 30 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 5 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 5 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 5 * ldG), gik);
      }
   }

   // cleanup code
   for(; p_outer < npts; p_outer += NPTS_LOCAL) {
      size_t npts_inner = MIN((size_t) NPTS_LOCAL, npts - p_outer);
      double *_point_outer = (_points + p_outer);

      double xA = shpair.rA.x;
      double yA = shpair.rA.y;
      double zA = shpair.rA.z;

      for(int i = 0; i < 31 * NPTS_LOCAL; i += SIMD_LENGTH) SIMD_ALIGNED_STORE((temp + i), SIMD_ZERO());

      for(int ij = 0; ij < shpair.nprim_pair; ++ij) {
         double RHO = shpair.prim_pairs[ij].gamma;
         double RHO_INV = 1.0 / RHO;

         constexpr double X_PA = 0.0;
         constexpr double Y_PA = 0.0;
         constexpr double Z_PA = 0.0;

         double eval = shpair.prim_pairs[ij].coeff_prod * shpair.prim_pairs[ij].K;

         // Evaluate T Values
         size_t npts_inner_upper = SIMD_LENGTH * (npts_inner / SIMD_LENGTH);
         size_t p_inner = 0;
         for(p_inner = 0; p_inner < npts_inner_upper; p_inner += SIMD_LENGTH) {
            SIMD_TYPE xC = SIMD_UNALIGNED_LOAD((_point_outer + p_inner + 0 * npts));
            SIMD_TYPE yC = SIMD_UNALIGNED_LOAD((_point_outer + p_inner + 1 * npts));
            SIMD_TYPE zC = SIMD_UNALIGNED_LOAD((_point_outer + p_inner + 2 * npts));

            SIMD_TYPE X_PC = SIMD_SUB(SIMD_DUPLICATE(&(xA)), xC);
            SIMD_TYPE Y_PC = SIMD_SUB(SIMD_DUPLICATE(&(yA)), yC);
            SIMD_TYPE Z_PC = SIMD_SUB(SIMD_DUPLICATE(&(zA)), zC);

            X_PC = SIMD_MUL(X_PC, X_PC);
            X_PC = SIMD_FMA(Y_PC, Y_PC, X_PC);
            X_PC = SIMD_FMA(Z_PC, Z_PC, X_PC);
            X_PC = SIMD_MUL(SIMD_DUPLICATE(&(RHO)), X_PC);
            SIMD_ALIGNED_STORE((Tval + p_inner), X_PC);
         }

         for(; p_inner < npts_inner; p_inner += SCALAR_LENGTH) {
            SCALAR_TYPE xC = SCALAR_LOAD((_point_outer + p_inner + 0 * npts));
            SCALAR_TYPE yC = SCALAR_LOAD((_point_outer + p_inner + 1 * npts));
            SCALAR_TYPE zC = SCALAR_LOAD((_point_outer + p_inner + 2 * npts));

            SCALAR_TYPE X_PC = SCALAR_SUB(SCALAR_DUPLICATE(&(xA)), xC);
            SCALAR_TYPE Y_PC = SCALAR_SUB(SCALAR_DUPLICATE(&(yA)), yC);
            SCALAR_TYPE Z_PC = SCALAR_SUB(SCALAR_DUPLICATE(&(zA)), zC);

            X_PC = SCALAR_MUL(X_PC, X_PC);
            X_PC = SCALAR_FMA(Y_PC, Y_PC, X_PC);
            X_PC = SCALAR_FMA(Z_PC, Z_PC, X_PC);
            X_PC = SCALAR_MUL(SCALAR_DUPLICATE(&(RHO)), X_PC);
            SCALAR_STORE((Tval + p_inner), X_PC);
         }

         // Evaluate Boys function
         GauXC::boys_function<0>(npts_inner, Tval, FmT + 0 * NPTS_LOCAL);
         GauXC::boys_function<1>(npts_inner, Tval, FmT + 1 * NPTS_LOCAL);
         GauXC::boys_function<2>(npts_inner, Tval, FmT + 2 * NPTS_LOCAL);
         GauXC::boys_function<3>(npts_inner, Tval, FmT + 3 * NPTS_LOCAL);
         GauXC::boys_function<4>(npts_inner, Tval, FmT + 4 * NPTS_LOCAL);

         // Evaluate VRR Buffer
         p_inner = 0;
         for(p_inner = 0; p_inner < npts_inner_upper; p_inner += SIMD_LENGTH) {
            SIMD_TYPE xC = SIMD_UNALIGNED_LOAD((_point_outer + p_inner + 0 * npts));
            SIMD_TYPE yC = SIMD_UNALIGNED_LOAD((_point_outer + p_inner + 1 * npts));
            SIMD_TYPE zC = SIMD_UNALIGNED_LOAD((_point_outer + p_inner + 2 * npts));

            SIMD_TYPE X_PC = SIMD_SUB(SIMD_DUPLICATE(&(xA)), xC);
            SIMD_TYPE Y_PC = SIMD_SUB(SIMD_DUPLICATE(&(yA)), yC);
            SIMD_TYPE Z_PC = SIMD_SUB(SIMD_DUPLICATE(&(zA)), zC);

            SIMD_TYPE tx, ty, t00, t01, t02, t03, t04, t10, t11, t12, t13, t20, t21, t22, t30, t31, t40;

            t00 = SIMD_ALIGNED_LOAD((FmT + p_inner + 0 * NPTS_LOCAL));
            t00 = SIMD_MUL(SIMD_DUPLICATE(&(eval)), t00);
            t01 = SIMD_ALIGNED_LOAD((FmT + p_inner + 1 * NPTS_LOCAL));
            t01 = SIMD_MUL(SIMD_DUPLICATE(&(eval)), t01);
            t02 = SIMD_ALIGNED_LOAD((FmT + p_inner + 2 * NPTS_LOCAL));
            t02 = SIMD_MUL(SIMD_DUPLICATE(&(eval)), t02);
            t03 = SIMD_ALIGNED_LOAD((FmT + p_inner + 3 * NPTS_LOCAL));
            t03 = SIMD_MUL(SIMD_DUPLICATE(&(eval)), t03);
            t04 = SIMD_ALIGNED_LOAD((FmT + p_inner + 4 * NPTS_LOCAL));
            t04 = SIMD_MUL(SIMD_DUPLICATE(&(eval)), t04);
            t10 = SIMD_MUL(SIMD_DUPLICATE(&(X_PA)), t00);
            t10 = SIMD_FNMA(X_PC, t01, t10);
            t11 = SIMD_MUL(SIMD_DUPLICATE(&(X_PA)), t01);
            t11 = SIMD_FNMA(X_PC, t02, t11);
            t12 = SIMD_MUL(SIMD_DUPLICATE(&(X_PA)), t02);
            t12 = SIMD_FNMA(X_PC, t03, t12);
            t13 = SIMD_MUL(SIMD_DUPLICATE(&(X_PA)), t03);
            t13 = SIMD_FNMA(X_PC, t04, t13);
            t20 = SIMD_MUL(SIMD_DUPLICATE(&(X_PA)), t10);
            t20 = SIMD_FNMA(X_PC, t11, t20);
            tx = SIMD_SUB(t00, t01);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t20 = SIMD_FMA(tx, ty, t20);
            t21 = SIMD_MUL(SIMD_DUPLICATE(&(X_PA)), t11);
            t21 = SIMD_FNMA(X_PC, t12, t21);
            tx = SIMD_SUB(t01, t02);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t21 = SIMD_FMA(tx, ty, t21);
            t22 = SIMD_MUL(SIMD_DUPLICATE(&(X_PA)), t12);
            t22 = SIMD_FNMA(X_PC, t13, t22);
            tx = SIMD_SUB(t02, t03);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t22 = SIMD_FMA(tx, ty, t22);
            tx = SIMD_ALIGNED_LOAD((temp + 0 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t20);
            SIMD_ALIGNED_STORE((temp + 0 * NPTS_LOCAL + p_inner), tx);
            t30 = SIMD_MUL(SIMD_DUPLICATE(&(X_PA)), t20);
            t30 = SIMD_FNMA(X_PC, t21, t30);
            tx = SIMD_SUB(t10, t11);
            ty = SIMD_SET1(0.5 * 2);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t30 = SIMD_FMA(tx, ty, t30);
            t31 = SIMD_MUL(SIMD_DUPLICATE(&(X_PA)), t21);
            t31 = SIMD_FNMA(X_PC, t22, t31);
            tx = SIMD_SUB(t11, t12);
            ty = SIMD_SET1(0.5 * 2);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t31 = SIMD_FMA(tx, ty, t31);
            tx = SIMD_ALIGNED_LOAD((temp + 6 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t30);
            SIMD_ALIGNED_STORE((temp + 6 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(X_PA)), t30);
            t40 = SIMD_FNMA(X_PC, t31, t40);
            tx = SIMD_SUB(t20, t21);
            ty = SIMD_SET1(0.5 * 3);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t40 = SIMD_FMA(tx, ty, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 16 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 16 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t30);
            t40 = SIMD_FNMA(Y_PC, t31, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 17 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 17 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t30);
            t40 = SIMD_FNMA(Z_PC, t31, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 18 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 18 * NPTS_LOCAL + p_inner), tx);
            t30 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t20);
            t30 = SIMD_FNMA(Y_PC, t21, t30);
            t31 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t21);
            t31 = SIMD_FNMA(Y_PC, t22, t31);
            tx = SIMD_ALIGNED_LOAD((temp + 7 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t30);
            SIMD_ALIGNED_STORE((temp + 7 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t30);
            t40 = SIMD_FNMA(Y_PC, t31, t40);
            tx = SIMD_SUB(t20, t21);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t40 = SIMD_FMA(tx, ty, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 19 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 19 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t30);
            t40 = SIMD_FNMA(Z_PC, t31, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 20 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 20 * NPTS_LOCAL + p_inner), tx);
            t30 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t20);
            t30 = SIMD_FNMA(Z_PC, t21, t30);
            t31 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t21);
            t31 = SIMD_FNMA(Z_PC, t22, t31);
            tx = SIMD_ALIGNED_LOAD((temp + 8 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t30);
            SIMD_ALIGNED_STORE((temp + 8 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t30);
            t40 = SIMD_FNMA(Z_PC, t31, t40);
            tx = SIMD_SUB(t20, t21);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t40 = SIMD_FMA(tx, ty, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 21 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 21 * NPTS_LOCAL + p_inner), tx);
            t20 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t10);
            t20 = SIMD_FNMA(Y_PC, t11, t20);
            t21 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t11);
            t21 = SIMD_FNMA(Y_PC, t12, t21);
            t22 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t12);
            t22 = SIMD_FNMA(Y_PC, t13, t22);
            tx = SIMD_ALIGNED_LOAD((temp + 1 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t20);
            SIMD_ALIGNED_STORE((temp + 1 * NPTS_LOCAL + p_inner), tx);
            t30 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t20);
            t30 = SIMD_FNMA(Y_PC, t21, t30);
            tx = SIMD_SUB(t10, t11);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t30 = SIMD_FMA(tx, ty, t30);
            t31 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t21);
            t31 = SIMD_FNMA(Y_PC, t22, t31);
            tx = SIMD_SUB(t11, t12);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t31 = SIMD_FMA(tx, ty, t31);
            tx = SIMD_ALIGNED_LOAD((temp + 9 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t30);
            SIMD_ALIGNED_STORE((temp + 9 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t30);
            t40 = SIMD_FNMA(Y_PC, t31, t40);
            tx = SIMD_SUB(t20, t21);
            ty = SIMD_SET1(0.5 * 2);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t40 = SIMD_FMA(tx, ty, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 22 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 22 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t30);
            t40 = SIMD_FNMA(Z_PC, t31, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 23 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 23 * NPTS_LOCAL + p_inner), tx);
            t30 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t20);
            t30 = SIMD_FNMA(Z_PC, t21, t30);
            t31 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t21);
            t31 = SIMD_FNMA(Z_PC, t22, t31);
            tx = SIMD_ALIGNED_LOAD((temp + 10 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t30);
            SIMD_ALIGNED_STORE((temp + 10 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t30);
            t40 = SIMD_FNMA(Z_PC, t31, t40);
            tx = SIMD_SUB(t20, t21);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t40 = SIMD_FMA(tx, ty, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 24 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 24 * NPTS_LOCAL + p_inner), tx);
            t20 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t10);
            t20 = SIMD_FNMA(Z_PC, t11, t20);
            t21 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t11);
            t21 = SIMD_FNMA(Z_PC, t12, t21);
            t22 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t12);
            t22 = SIMD_FNMA(Z_PC, t13, t22);
            tx = SIMD_ALIGNED_LOAD((temp + 2 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t20);
            SIMD_ALIGNED_STORE((temp + 2 * NPTS_LOCAL + p_inner), tx);
            t30 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t20);
            t30 = SIMD_FNMA(Z_PC, t21, t30);
            tx = SIMD_SUB(t10, t11);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t30 = SIMD_FMA(tx, ty, t30);
            t31 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t21);
            t31 = SIMD_FNMA(Z_PC, t22, t31);
            tx = SIMD_SUB(t11, t12);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t31 = SIMD_FMA(tx, ty, t31);
            tx = SIMD_ALIGNED_LOAD((temp + 11 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t30);
            SIMD_ALIGNED_STORE((temp + 11 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t30);
            t40 = SIMD_FNMA(Z_PC, t31, t40);
            tx = SIMD_SUB(t20, t21);
            ty = SIMD_SET1(0.5 * 2);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t40 = SIMD_FMA(tx, ty, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 25 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 25 * NPTS_LOCAL + p_inner), tx);
            t10 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t00);
            t10 = SIMD_FNMA(Y_PC, t01, t10);
            t11 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t01);
            t11 = SIMD_FNMA(Y_PC, t02, t11);
            t12 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t02);
            t12 = SIMD_FNMA(Y_PC, t03, t12);
            t13 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t03);
            t13 = SIMD_FNMA(Y_PC, t04, t13);
            t20 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t10);
            t20 = SIMD_FNMA(Y_PC, t11, t20);
            tx = SIMD_SUB(t00, t01);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t20 = SIMD_FMA(tx, ty, t20);
            t21 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t11);
            t21 = SIMD_FNMA(Y_PC, t12, t21);
            tx = SIMD_SUB(t01, t02);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t21 = SIMD_FMA(tx, ty, t21);
            t22 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t12);
            t22 = SIMD_FNMA(Y_PC, t13, t22);
            tx = SIMD_SUB(t02, t03);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t22 = SIMD_FMA(tx, ty, t22);
            tx = SIMD_ALIGNED_LOAD((temp + 3 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t20);
            SIMD_ALIGNED_STORE((temp + 3 * NPTS_LOCAL + p_inner), tx);
            t30 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t20);
            t30 = SIMD_FNMA(Y_PC, t21, t30);
            tx = SIMD_SUB(t10, t11);
            ty = SIMD_SET1(0.5 * 2);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t30 = SIMD_FMA(tx, ty, t30);
            t31 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t21);
            t31 = SIMD_FNMA(Y_PC, t22, t31);
            tx = SIMD_SUB(t11, t12);
            ty = SIMD_SET1(0.5 * 2);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t31 = SIMD_FMA(tx, ty, t31);
            tx = SIMD_ALIGNED_LOAD((temp + 12 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t30);
            SIMD_ALIGNED_STORE((temp + 12 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(Y_PA)), t30);
            t40 = SIMD_FNMA(Y_PC, t31, t40);
            tx = SIMD_SUB(t20, t21);
            ty = SIMD_SET1(0.5 * 3);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t40 = SIMD_FMA(tx, ty, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 26 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 26 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t30);
            t40 = SIMD_FNMA(Z_PC, t31, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 27 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 27 * NPTS_LOCAL + p_inner), tx);
            t30 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t20);
            t30 = SIMD_FNMA(Z_PC, t21, t30);
            t31 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t21);
            t31 = SIMD_FNMA(Z_PC, t22, t31);
            tx = SIMD_ALIGNED_LOAD((temp + 13 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t30);
            SIMD_ALIGNED_STORE((temp + 13 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t30);
            t40 = SIMD_FNMA(Z_PC, t31, t40);
            tx = SIMD_SUB(t20, t21);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t40 = SIMD_FMA(tx, ty, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 28 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 28 * NPTS_LOCAL + p_inner), tx);
            t20 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t10);
            t20 = SIMD_FNMA(Z_PC, t11, t20);
            t21 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t11);
            t21 = SIMD_FNMA(Z_PC, t12, t21);
            t22 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t12);
            t22 = SIMD_FNMA(Z_PC, t13, t22);
            tx = SIMD_ALIGNED_LOAD((temp + 4 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t20);
            SIMD_ALIGNED_STORE((temp + 4 * NPTS_LOCAL + p_inner), tx);
            t30 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t20);
            t30 = SIMD_FNMA(Z_PC, t21, t30);
            tx = SIMD_SUB(t10, t11);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t30 = SIMD_FMA(tx, ty, t30);
            t31 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t21);
            t31 = SIMD_FNMA(Z_PC, t22, t31);
            tx = SIMD_SUB(t11, t12);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t31 = SIMD_FMA(tx, ty, t31);
            tx = SIMD_ALIGNED_LOAD((temp + 14 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t30);
            SIMD_ALIGNED_STORE((temp + 14 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t30);
            t40 = SIMD_FNMA(Z_PC, t31, t40);
            tx = SIMD_SUB(t20, t21);
            ty = SIMD_SET1(0.5 * 2);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t40 = SIMD_FMA(tx, ty, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 29 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 29 * NPTS_LOCAL + p_inner), tx);
            t10 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t00);
            t10 = SIMD_FNMA(Z_PC, t01, t10);
            t11 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t01);
            t11 = SIMD_FNMA(Z_PC, t02, t11);
            t12 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t02);
            t12 = SIMD_FNMA(Z_PC, t03, t12);
            t13 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t03);
            t13 = SIMD_FNMA(Z_PC, t04, t13);
            t20 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t10);
            t20 = SIMD_FNMA(Z_PC, t11, t20);
            tx = SIMD_SUB(t00, t01);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t20 = SIMD_FMA(tx, ty, t20);
            t21 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t11);
            t21 = SIMD_FNMA(Z_PC, t12, t21);
            tx = SIMD_SUB(t01, t02);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t21 = SIMD_FMA(tx, ty, t21);
            t22 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t12);
            t22 = SIMD_FNMA(Z_PC, t13, t22);
            tx = SIMD_SUB(t02, t03);
            ty = SIMD_SET1(0.5 * 1);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t22 = SIMD_FMA(tx, ty, t22);
            tx = SIMD_ALIGNED_LOAD((temp + 5 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t20);
            SIMD_ALIGNED_STORE((temp + 5 * NPTS_LOCAL + p_inner), tx);
            t30 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t20);
            t30 = SIMD_FNMA(Z_PC, t21, t30);
            tx = SIMD_SUB(t10, t11);
            ty = SIMD_SET1(0.5 * 2);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t30 = SIMD_FMA(tx, ty, t30);
            t31 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t21);
            t31 = SIMD_FNMA(Z_PC, t22, t31);
            tx = SIMD_SUB(t11, t12);
            ty = SIMD_SET1(0.5 * 2);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t31 = SIMD_FMA(tx, ty, t31);
            tx = SIMD_ALIGNED_LOAD((temp + 15 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t30);
            SIMD_ALIGNED_STORE((temp + 15 * NPTS_LOCAL + p_inner), tx);
            t40 = SIMD_MUL(SIMD_DUPLICATE(&(Z_PA)), t30);
            t40 = SIMD_FNMA(Z_PC, t31, t40);
            tx = SIMD_SUB(t20, t21);
            ty = SIMD_SET1(0.5 * 3);
            ty = SIMD_MUL(ty, SIMD_DUPLICATE(&(RHO_INV)));
            t40 = SIMD_FMA(tx, ty, t40);
            tx = SIMD_ALIGNED_LOAD((temp + 30 * NPTS_LOCAL + p_inner));
            tx = SIMD_ADD(tx, t40);
            SIMD_ALIGNED_STORE((temp + 30 * NPTS_LOCAL + p_inner), tx);
         }

         for(; p_inner < npts_inner; p_inner += SCALAR_LENGTH) {
            SCALAR_TYPE xC = SCALAR_LOAD((_point_outer + p_inner + 0 * npts));
            SCALAR_TYPE yC = SCALAR_LOAD((_point_outer + p_inner + 1 * npts));
            SCALAR_TYPE zC = SCALAR_LOAD((_point_outer + p_inner + 2 * npts));

            SCALAR_TYPE X_PC = SCALAR_SUB(SCALAR_DUPLICATE(&(xA)), xC);
            SCALAR_TYPE Y_PC = SCALAR_SUB(SCALAR_DUPLICATE(&(yA)), yC);
            SCALAR_TYPE Z_PC = SCALAR_SUB(SCALAR_DUPLICATE(&(zA)), zC);

            SCALAR_TYPE tx, ty, t00, t01, t02, t03, t04, t10, t11, t12, t13, t20, t21, t22, t30, t31, t40;

            t00 = SCALAR_LOAD((FmT + p_inner + 0 * NPTS_LOCAL));
            t00 = SCALAR_MUL(SCALAR_DUPLICATE(&(eval)), t00);
            t01 = SCALAR_LOAD((FmT + p_inner + 1 * NPTS_LOCAL));
            t01 = SCALAR_MUL(SCALAR_DUPLICATE(&(eval)), t01);
            t02 = SCALAR_LOAD((FmT + p_inner + 2 * NPTS_LOCAL));
            t02 = SCALAR_MUL(SCALAR_DUPLICATE(&(eval)), t02);
            t03 = SCALAR_LOAD((FmT + p_inner + 3 * NPTS_LOCAL));
            t03 = SCALAR_MUL(SCALAR_DUPLICATE(&(eval)), t03);
            t04 = SCALAR_LOAD((FmT + p_inner + 4 * NPTS_LOCAL));
            t04 = SCALAR_MUL(SCALAR_DUPLICATE(&(eval)), t04);
            t10 = SCALAR_MUL(SCALAR_DUPLICATE(&(X_PA)), t00);
            t10 = SCALAR_FNMA(X_PC, t01, t10);
            t11 = SCALAR_MUL(SCALAR_DUPLICATE(&(X_PA)), t01);
            t11 = SCALAR_FNMA(X_PC, t02, t11);
            t12 = SCALAR_MUL(SCALAR_DUPLICATE(&(X_PA)), t02);
            t12 = SCALAR_FNMA(X_PC, t03, t12);
            t13 = SCALAR_MUL(SCALAR_DUPLICATE(&(X_PA)), t03);
            t13 = SCALAR_FNMA(X_PC, t04, t13);
            t20 = SCALAR_MUL(SCALAR_DUPLICATE(&(X_PA)), t10);
            t20 = SCALAR_FNMA(X_PC, t11, t20);
            tx = SCALAR_SUB(t00, t01);
            ty = SCALAR_SET1(0.5 * 1);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t20 = SCALAR_FMA(tx, ty, t20);
            t21 = SCALAR_MUL(SCALAR_DUPLICATE(&(X_PA)), t11);
            t21 = SCALAR_FNMA(X_PC, t12, t21);
            tx = SCALAR_SUB(t01, t02);
            ty = SCALAR_SET1(0.5 * 1);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t21 = SCALAR_FMA(tx, ty, t21);
            t22 = SCALAR_MUL(SCALAR_DUPLICATE(&(X_PA)), t12);
            t22 = SCALAR_FNMA(X_PC, t13, t22);
            tx = SCALAR_SUB(t02, t03);
            ty = SCALAR_SET1(0.5 * 1);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t22 = SCALAR_FMA(tx, ty, t22);
            tx = SCALAR_LOAD((temp + 0 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t20);
            SCALAR_STORE((temp + 0 * NPTS_LOCAL + p_inner), tx);
            t30 = SCALAR_MUL(SCALAR_DUPLICATE(&(X_PA)), t20);
            t30 = SCALAR_FNMA(X_PC, t21, t30);
            tx = SCALAR_SUB(t10, t11);
            ty = SCALAR_SET1(0.5 * 2);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t30 = SCALAR_FMA(tx, ty, t30);
            t31 = SCALAR_MUL(SCALAR_DUPLICATE(&(X_PA)), t21);
            t31 = SCALAR_FNMA(X_PC, t22, t31);
            tx = SCALAR_SUB(t11, t12);
            ty = SCALAR_SET1(0.5 * 2);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t31 = SCALAR_FMA(tx, ty, t31);
            tx = SCALAR_LOAD((temp + 6 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t30);
            SCALAR_STORE((temp + 6 * NPTS_LOCAL + p_inner), tx);
            t40 = SCALAR_MUL(SCALAR_DUPLICATE(&(X_PA)), t30);
            t40 = SCALAR_FNMA(X_PC, t31, t40);
            tx = SCALAR_SUB(t20, t21);
            ty = SCALAR_SET1(0.5 * 3);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t40 = SCALAR_FMA(tx, ty, t40);
            tx = SCALAR_LOAD((temp + 16 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t40);
            SCALAR_STORE((temp + 16 * NPTS_LOCAL + p_inner), tx);
            t40 = SCALAR_MUL(SCALAR_DUPLICATE(&(Y_PA)), t30);
            t40 = SCALAR_FNMA(Y_PC, t31, t40);
            tx = SCALAR_LOAD((temp + 17 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t40);
            SCALAR_STORE((temp + 17 * NPTS_LOCAL + p_inner), tx);
            t40 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t30);
            t40 = SCALAR_FNMA(Z_PC, t31, t40);
            tx = SCALAR_LOAD((temp + 18 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t40);
            SCALAR_STORE((temp + 18 * NPTS_LOCAL + p_inner), tx);
            t30 = SCALAR_MUL(SCALAR_DUPLICATE(&(Y_PA)), t20);
            t30 = SCALAR_FNMA(Y_PC, t21, t30);
            t31 = SCALAR_MUL(SCALAR_DUPLICATE(&(Y_PA)), t21);
            t31 = SCALAR_FNMA(Y_PC, t22, t31);
            tx = SCALAR_LOAD((temp + 7 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t30);
            SCALAR_STORE((temp + 7 * NPTS_LOCAL + p_inner), tx);
            t40 = SCALAR_MUL(SCALAR_DUPLICATE(&(Y_PA)), t30);
            t40 = SCALAR_FNMA(Y_PC, t31, t40);
            tx = SCALAR_SUB(t20, t21);
            ty = SCALAR_SET1(0.5 * 1);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t40 = SCALAR_FMA(tx, ty, t40);
            tx = SCALAR_LOAD((temp + 19 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t40);
            SCALAR_STORE((temp + 19 * NPTS_LOCAL + p_inner), tx);
            t40 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t30);
            t40 = SCALAR_FNMA(Z_PC, t31, t40);
            tx = SCALAR_LOAD((temp + 20 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t40);
            SCALAR_STORE((temp + 20 * NPTS_LOCAL + p_inner), tx);
            t30 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t20);
            t30 = SCALAR_FNMA(Z_PC, t21, t30);
            t31 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t21);
            t31 = SCALAR_FNMA(Z_PC, t22, t31);
            tx = SCALAR_LOAD((temp + 8 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t30);
            SCALAR_STORE((temp + 8 * NPTS_LOCAL + p_inner), tx);
            t40 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t30);
            t40 = SCALAR_FNMA(Z_PC, t31, t40);
            tx = SCALAR_SUB(t20, t21);
            ty = SCALAR_SET1(0.5 * 1);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t40 = SCALAR_FMA(tx, ty, t40);
            tx = SCALAR_LOAD((temp + 21 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t40);
            SCALAR_STORE((temp + 21 * NPTS_LOCAL + p_inner), tx);
            t20 = SCALAR_MUL(SCALAR_DUPLICATE(&(Y_PA)), t10);
            t20 = SCALAR_FNMA(Y_PC, t11, t20);
            t21 = SCALAR_MUL(SCALAR_DUPLICATE(&(Y_PA)), t11);
            t21 = SCALAR_FNMA(Y_PC, t12, t21);
            t22 = SCALAR_MUL(SCALAR_DUPLICATE(&(Y_PA)), t12);
            t22 = SCALAR_FNMA(Y_PC, t13, t22);
            tx = SCALAR_LOAD((temp + 1 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t20);
            SCALAR_STORE((temp + 1 * NPTS_LOCAL + p_inner), tx);
            t30 = SCALAR_MUL(SCALAR_DUPLICATE(&(Y_PA)), t20);
            t30 = SCALAR_FNMA(Y_PC, t21, t30);
            tx = SCALAR_SUB(t10, t11);
            ty = SCALAR_SET1(0.5 * 1);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t30 = SCALAR_FMA(tx, ty, t30);
            t31 = SCALAR_MUL(SCALAR_DUPLICATE(&(Y_PA)), t21);
            t31 = SCALAR_FNMA(Y_PC, t22, t31);
            tx = SCALAR_SUB(t11, t12);
            ty = SCALAR_SET1(0.5 * 1);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t31 = SCALAR_FMA(tx, ty, t31);
            tx = SCALAR_LOAD((temp + 9 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t30);
            SCALAR_STORE((temp + 9 * NPTS_LOCAL + p_inner), tx);
            t40 = SCALAR_MUL(SCALAR_DUPLICATE(&(Y_PA)), t30);
            t40 = SCALAR_FNMA(Y_PC, t31, t40);
            tx = SCALAR_SUB(t20, t21);
            ty = SCALAR_SET1(0.5 * 2);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t40 = SCALAR_FMA(tx, ty, t40);
            tx = SCALAR_LOAD((temp + 22 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t40);
            SCALAR_STORE((temp + 22 * NPTS_LOCAL + p_inner), tx);
            t40 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t30);
            t40 = SCALAR_FNMA(Z_PC, t31, t40);
            tx = SCALAR_LOAD((temp + 23 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t40);
            SCALAR_STORE((temp + 23 * NPTS_LOCAL + p_inner), tx);
            t30 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t20);
            t30 = SCALAR_FNMA(Z_PC, t21, t30);
            t31 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t21);
            t31 = SCALAR_FNMA(Z_PC, t22, t31);
            tx = SCALAR_LOAD((temp + 10 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t30);
            SCALAR_STORE((temp + 10 * NPTS_LOCAL + p_inner), tx);
            t40 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t30);
            t40 = SCALAR_FNMA(Z_PC, t31, t40);
            tx = SCALAR_SUB(t20, t21);
            ty = SCALAR_SET1(0.5 * 1);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t40 = SCALAR_FMA(tx, ty, t40);
            tx = SCALAR_LOAD((temp + 24 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t40);
            SCALAR_STORE((temp + 24 * NPTS_LOCAL + p_inner), tx);
            t20 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t10);
            t20 = SCALAR_FNMA(Z_PC, t11, t20);
            t21 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t11);
            t21 = SCALAR_FNMA(Z_PC, t12, t21);
            t22 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t12);
            t22 = SCALAR_FNMA(Z_PC, t13, t22);
            tx = SCALAR_LOAD((temp + 2 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t20);
            SCALAR_STORE((temp + 2 * NPTS_LOCAL + p_inner), tx);
            t30 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t20);
            t30 = SCALAR_FNMA(Z_PC, t21, t30);
            tx = SCALAR_SUB(t10, t11);
            ty = SCALAR_SET1(0.5 * 1);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t30 = SCALAR_FMA(tx, ty, t30);
            t31 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t21);
            t31 = SCALAR_FNMA(Z_PC, t22, t31);
            tx = SCALAR_SUB(t11, t12);
            ty = SCALAR_SET1(0.5 * 1);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t31 = SCALAR_FMA(tx, ty, t31);
            tx = SCALAR_LOAD((temp + 11 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t30);
            SCALAR_STORE((temp + 11 * NPTS_LOCAL + p_inner), tx);
            t40 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t30);
            t40 = SCALAR_FNMA(Z_PC, t31, t40);
            tx = SCALAR_SUB(t20, t21);
            ty = SCALAR_SET1(0.5 * 2);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t40 = SCALAR_FMA(tx, ty, t40);
            tx = SCALAR_LOAD((temp + 25 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t40);
            SCALAR_STORE((temp + 25 * NPTS_LOCAL + p_inner), tx);
            t10 = SCALAR_MUL(SCALAR_DUPLICATE(&(Y_PA)), t00);
            t10 = SCALAR_FNMA(Y_PC, t01, t10);
            t11 = SCALAR_MUL(SCALAR_DUPLICATE(&(Y_PA)), t01);
            t11 = SCALAR_FNMA(Y_PC, t02, t11);
            t12 = SCALAR_MUL(SCALAR_DUPLICATE(&(Y_PA)), t02);
            t12 = SCALAR_FNMA(Y_PC, t03, t12);
            t13 = SCALAR_MUL(SCALAR_DUPLICATE(&(Y_PA)), t03);
            t13 = SCALAR_FNMA(Y_PC, t04, t13);
            t20 = SCALAR_MUL(SCALAR_DUPLICATE(&(Y_PA)), t10);
            t20 = SCALAR_FNMA(Y_PC, t11, t20);
            tx = SCALAR_SUB(t00, t01);
            ty = SCALAR_SET1(0.5 * 1);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t20 = SCALAR_FMA(tx, ty, t20);
            t21 = SCALAR_MUL(SCALAR_DUPLICATE(&(Y_PA)), t11);
            t21 = SCALAR_FNMA(Y_PC, t12, t21);
            tx = SCALAR_SUB(t01, t02);
            ty = SCALAR_SET1(0.5 * 1);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t21 = SCALAR_FMA(tx, ty, t21);
            t22 = SCALAR_MUL(SCALAR_DUPLICATE(&(Y_PA)), t12);
            t22 = SCALAR_FNMA(Y_PC, t13, t22);
            tx = SCALAR_SUB(t02, t03);
            ty = SCALAR_SET1(0.5 * 1);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t22 = SCALAR_FMA(tx, ty, t22);
            tx = SCALAR_LOAD((temp + 3 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t20);
            SCALAR_STORE((temp + 3 * NPTS_LOCAL + p_inner), tx);
            t30 = SCALAR_MUL(SCALAR_DUPLICATE(&(Y_PA)), t20);
            t30 = SCALAR_FNMA(Y_PC, t21, t30);
            tx = SCALAR_SUB(t10, t11);
            ty = SCALAR_SET1(0.5 * 2);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t30 = SCALAR_FMA(tx, ty, t30);
            t31 = SCALAR_MUL(SCALAR_DUPLICATE(&(Y_PA)), t21);
            t31 = SCALAR_FNMA(Y_PC, t22, t31);
            tx = SCALAR_SUB(t11, t12);
            ty = SCALAR_SET1(0.5 * 2);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t31 = SCALAR_FMA(tx, ty, t31);
            tx = SCALAR_LOAD((temp + 12 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t30);
            SCALAR_STORE((temp + 12 * NPTS_LOCAL + p_inner), tx);
            t40 = SCALAR_MUL(SCALAR_DUPLICATE(&(Y_PA)), t30);
            t40 = SCALAR_FNMA(Y_PC, t31, t40);
            tx = SCALAR_SUB(t20, t21);
            ty = SCALAR_SET1(0.5 * 3);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t40 = SCALAR_FMA(tx, ty, t40);
            tx = SCALAR_LOAD((temp + 26 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t40);
            SCALAR_STORE((temp + 26 * NPTS_LOCAL + p_inner), tx);
            t40 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t30);
            t40 = SCALAR_FNMA(Z_PC, t31, t40);
            tx = SCALAR_LOAD((temp + 27 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t40);
            SCALAR_STORE((temp + 27 * NPTS_LOCAL + p_inner), tx);
            t30 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t20);
            t30 = SCALAR_FNMA(Z_PC, t21, t30);
            t31 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t21);
            t31 = SCALAR_FNMA(Z_PC, t22, t31);
            tx = SCALAR_LOAD((temp + 13 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t30);
            SCALAR_STORE((temp + 13 * NPTS_LOCAL + p_inner), tx);
            t40 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t30);
            t40 = SCALAR_FNMA(Z_PC, t31, t40);
            tx = SCALAR_SUB(t20, t21);
            ty = SCALAR_SET1(0.5 * 1);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t40 = SCALAR_FMA(tx, ty, t40);
            tx = SCALAR_LOAD((temp + 28 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t40);
            SCALAR_STORE((temp + 28 * NPTS_LOCAL + p_inner), tx);
            t20 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t10);
            t20 = SCALAR_FNMA(Z_PC, t11, t20);
            t21 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t11);
            t21 = SCALAR_FNMA(Z_PC, t12, t21);
            t22 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t12);
            t22 = SCALAR_FNMA(Z_PC, t13, t22);
            tx = SCALAR_LOAD((temp + 4 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t20);
            SCALAR_STORE((temp + 4 * NPTS_LOCAL + p_inner), tx);
            t30 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t20);
            t30 = SCALAR_FNMA(Z_PC, t21, t30);
            tx = SCALAR_SUB(t10, t11);
            ty = SCALAR_SET1(0.5 * 1);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t30 = SCALAR_FMA(tx, ty, t30);
            t31 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t21);
            t31 = SCALAR_FNMA(Z_PC, t22, t31);
            tx = SCALAR_SUB(t11, t12);
            ty = SCALAR_SET1(0.5 * 1);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t31 = SCALAR_FMA(tx, ty, t31);
            tx = SCALAR_LOAD((temp + 14 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t30);
            SCALAR_STORE((temp + 14 * NPTS_LOCAL + p_inner), tx);
            t40 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t30);
            t40 = SCALAR_FNMA(Z_PC, t31, t40);
            tx = SCALAR_SUB(t20, t21);
            ty = SCALAR_SET1(0.5 * 2);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t40 = SCALAR_FMA(tx, ty, t40);
            tx = SCALAR_LOAD((temp + 29 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t40);
            SCALAR_STORE((temp + 29 * NPTS_LOCAL + p_inner), tx);
            t10 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t00);
            t10 = SCALAR_FNMA(Z_PC, t01, t10);
            t11 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t01);
            t11 = SCALAR_FNMA(Z_PC, t02, t11);
            t12 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t02);
            t12 = SCALAR_FNMA(Z_PC, t03, t12);
            t13 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t03);
            t13 = SCALAR_FNMA(Z_PC, t04, t13);
            t20 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t10);
            t20 = SCALAR_FNMA(Z_PC, t11, t20);
            tx = SCALAR_SUB(t00, t01);
            ty = SCALAR_SET1(0.5 * 1);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t20 = SCALAR_FMA(tx, ty, t20);
            t21 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t11);
            t21 = SCALAR_FNMA(Z_PC, t12, t21);
            tx = SCALAR_SUB(t01, t02);
            ty = SCALAR_SET1(0.5 * 1);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t21 = SCALAR_FMA(tx, ty, t21);
            t22 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t12);
            t22 = SCALAR_FNMA(Z_PC, t13, t22);
            tx = SCALAR_SUB(t02, t03);
            ty = SCALAR_SET1(0.5 * 1);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t22 = SCALAR_FMA(tx, ty, t22);
            tx = SCALAR_LOAD((temp + 5 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t20);
            SCALAR_STORE((temp + 5 * NPTS_LOCAL + p_inner), tx);
            t30 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t20);
            t30 = SCALAR_FNMA(Z_PC, t21, t30);
            tx = SCALAR_SUB(t10, t11);
            ty = SCALAR_SET1(0.5 * 2);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t30 = SCALAR_FMA(tx, ty, t30);
            t31 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t21);
            t31 = SCALAR_FNMA(Z_PC, t22, t31);
            tx = SCALAR_SUB(t11, t12);
            ty = SCALAR_SET1(0.5 * 2);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t31 = SCALAR_FMA(tx, ty, t31);
            tx = SCALAR_LOAD((temp + 15 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t30);
            SCALAR_STORE((temp + 15 * NPTS_LOCAL + p_inner), tx);
            t40 = SCALAR_MUL(SCALAR_DUPLICATE(&(Z_PA)), t30);
            t40 = SCALAR_FNMA(Z_PC, t31, t40);
            tx = SCALAR_SUB(t20, t21);
            ty = SCALAR_SET1(0.5 * 3);
            ty = SCALAR_MUL(ty, SCALAR_DUPLICATE(&(RHO_INV)));
            t40 = SCALAR_FMA(tx, ty, t40);
            tx = SCALAR_LOAD((temp + 30 * NPTS_LOCAL + p_inner));
            tx = SCALAR_ADD(tx, t40);
            SCALAR_STORE((temp + 30 * NPTS_LOCAL + p_inner), tx);
         }
      }

      size_t npts_inner_upper = SIMD_LENGTH * (npts_inner / SIMD_LENGTH);
      size_t p_inner = 0;
      for(p_inner = 0; p_inner < npts_inner_upper; p_inner += SIMD_LENGTH) {
         double *Xik = (Xi + p_outer + p_inner);
         double *Gik = (Gi + p_outer + p_inner);

         SIMD_TYPE tx, wg, xik, gik;
         tx  = SIMD_ALIGNED_LOAD((temp + 16 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 0 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 0 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 0 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 17 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 0 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 1 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 1 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 18 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 0 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 2 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 2 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 19 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 0 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 3 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 3 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 20 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 0 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 4 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 4 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 21 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 0 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 5 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 5 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 17 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 1 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 0 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 0 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 19 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 1 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 1 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 1 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 20 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 1 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 2 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 2 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 22 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 1 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 3 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 3 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 23 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 1 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 4 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 4 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 24 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 1 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 5 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 5 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 18 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 2 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 0 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 0 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 20 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 2 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 1 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 1 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 21 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 2 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 2 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 2 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 23 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 2 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 3 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 3 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 24 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 2 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 4 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 4 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 25 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 2 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 5 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 5 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 19 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 3 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 0 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 0 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 22 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 3 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 1 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 1 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 23 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 3 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 2 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 2 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 26 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 3 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 3 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 3 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 27 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 3 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 4 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 4 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 28 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 3 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 5 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 5 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 20 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 4 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 0 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 0 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 23 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 4 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 1 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 1 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 24 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 4 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 2 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 2 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 27 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 4 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 3 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 3 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 28 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 4 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 4 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 4 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 29 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 4 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 5 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 5 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 21 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 5 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 0 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 0 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 24 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 5 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 1 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 1 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 25 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 5 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 2 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 2 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 28 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 5 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 3 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 3 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 29 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 5 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 4 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 4 * ldG), gik);
         tx  = SIMD_ALIGNED_LOAD((temp + 30 * NPTS_LOCAL + p_inner));
         wg  = SIMD_UNALIGNED_LOAD((weights + p_outer + p_inner));

         xik = SIMD_UNALIGNED_LOAD((Xik + 5 * ldX));
         gik = SIMD_UNALIGNED_LOAD((Gik + 5 * ldG));

         tx = SIMD_MUL(tx, wg);
         gik = SIMD_FMA(tx, xik, gik);
         SIMD_UNALIGNED_STORE((Gik + 5 * ldG), gik);
      }

      for(; p_inner < npts_inner; p_inner += SCALAR_LENGTH) {
         double *Xik = (Xi + p_outer + p_inner);
         double *Gik = (Gi + p_outer + p_inner);

         SCALAR_TYPE tx, wg, xik, gik;
         tx  = SCALAR_LOAD((temp + 16 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 0 * ldX));
         gik = SCALAR_LOAD((Gik + 0 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 0 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 17 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 0 * ldX));
         gik = SCALAR_LOAD((Gik + 1 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 1 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 18 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 0 * ldX));
         gik = SCALAR_LOAD((Gik + 2 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 2 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 19 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 0 * ldX));
         gik = SCALAR_LOAD((Gik + 3 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 3 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 20 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 0 * ldX));
         gik = SCALAR_LOAD((Gik + 4 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 4 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 21 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 0 * ldX));
         gik = SCALAR_LOAD((Gik + 5 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 5 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 17 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 1 * ldX));
         gik = SCALAR_LOAD((Gik + 0 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 0 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 19 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 1 * ldX));
         gik = SCALAR_LOAD((Gik + 1 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 1 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 20 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 1 * ldX));
         gik = SCALAR_LOAD((Gik + 2 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 2 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 22 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 1 * ldX));
         gik = SCALAR_LOAD((Gik + 3 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 3 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 23 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 1 * ldX));
         gik = SCALAR_LOAD((Gik + 4 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 4 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 24 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 1 * ldX));
         gik = SCALAR_LOAD((Gik + 5 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 5 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 18 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 2 * ldX));
         gik = SCALAR_LOAD((Gik + 0 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 0 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 20 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 2 * ldX));
         gik = SCALAR_LOAD((Gik + 1 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 1 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 21 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 2 * ldX));
         gik = SCALAR_LOAD((Gik + 2 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 2 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 23 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 2 * ldX));
         gik = SCALAR_LOAD((Gik + 3 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 3 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 24 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 2 * ldX));
         gik = SCALAR_LOAD((Gik + 4 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 4 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 25 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 2 * ldX));
         gik = SCALAR_LOAD((Gik + 5 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 5 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 19 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 3 * ldX));
         gik = SCALAR_LOAD((Gik + 0 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 0 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 22 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 3 * ldX));
         gik = SCALAR_LOAD((Gik + 1 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 1 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 23 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 3 * ldX));
         gik = SCALAR_LOAD((Gik + 2 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 2 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 26 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 3 * ldX));
         gik = SCALAR_LOAD((Gik + 3 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 3 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 27 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 3 * ldX));
         gik = SCALAR_LOAD((Gik + 4 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 4 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 28 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 3 * ldX));
         gik = SCALAR_LOAD((Gik + 5 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 5 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 20 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 4 * ldX));
         gik = SCALAR_LOAD((Gik + 0 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 0 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 23 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 4 * ldX));
         gik = SCALAR_LOAD((Gik + 1 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 1 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 24 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 4 * ldX));
         gik = SCALAR_LOAD((Gik + 2 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 2 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 27 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 4 * ldX));
         gik = SCALAR_LOAD((Gik + 3 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 3 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 28 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 4 * ldX));
         gik = SCALAR_LOAD((Gik + 4 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 4 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 29 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 4 * ldX));
         gik = SCALAR_LOAD((Gik + 5 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 5 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 21 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 5 * ldX));
         gik = SCALAR_LOAD((Gik + 0 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 0 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 24 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 5 * ldX));
         gik = SCALAR_LOAD((Gik + 1 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 1 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 25 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 5 * ldX));
         gik = SCALAR_LOAD((Gik + 2 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 2 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 28 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 5 * ldX));
         gik = SCALAR_LOAD((Gik + 3 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 3 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 29 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 5 * ldX));
         gik = SCALAR_LOAD((Gik + 4 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 4 * ldG), gik);
         tx  = SCALAR_LOAD((temp + 30 * NPTS_LOCAL + p_inner));
         wg  = SCALAR_LOAD((weights + p_outer + p_inner));

         xik = SCALAR_LOAD((Xik + 5 * ldX));
         gik = SCALAR_LOAD((Gik + 5 * ldG));

         tx = SCALAR_MUL(tx, wg);
         gik = SCALAR_FMA(tx, xik, gik);
         SCALAR_STORE((Gik + 5 * ldG), gik);
      }
   }
}
