#include <math.h>
#include "../include/gpu/chebyshev_boys_computation.hpp"
#include "../include/gpu/integral_data_types.hpp"
#include "config_obara_saika.hpp"
#include "integral_2_2.hu"

#define PI 3.14159265358979323846

#define MIN(a,b)			\
  ({ __typeof__ (a) _a = (a);	        \
  __typeof__ (b) _b = (b);		\
  _a < _b ? _a : _b; })

namespace XGPU {
__global__ void integral_2_2(size_t npts,
                             double *_points,
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
                             double *boys_table) {
   __shared__ double *temp;
   for(size_t p_outer = blockIdx.x * blockDim.x; p_outer < npts; p_outer += gridDim.x * blockDim.x) {
      double *_point_outer = (_points + p_outer);

      size_t p_inner = (threadIdx.x < (npts - p_outer)) ? threadIdx.x : (npts - p_outer);

      double X_AB = rA.x - rB.x;
      double Y_AB = rA.y - rB.y;
      double Z_AB = rA.z - rB.z;

      for(int i = 0; i < 31; ++i) SCALAR_STORE((temp + i * blockDim.x + threadIdx.x), SCALAR_ZERO());

      for(int ij = 0; ij < nprim_pairs; ++ij) {
         double RHO = prim_pairs[ij].gamma;
         double RHO_INV = prim_pairs[ij].gamma_inv;
         double X_PA = prim_pairs[ij].PA.x;
         double Y_PA = prim_pairs[ij].PA.y;
         double Z_PA = prim_pairs[ij].PA.z;

         double xP = prim_pairs[ij].P.x;
         double yP = prim_pairs[ij].P.y;
         double zP = prim_pairs[ij].P.z;

         double eval = prim_pairs[ij].K_coeff_prod;

         // Evaluate T Values
         SCALAR_TYPE xC = SCALAR_LOAD((_point_outer + p_inner + 0 * npts));
         SCALAR_TYPE yC = SCALAR_LOAD((_point_outer + p_inner + 1 * npts));
         SCALAR_TYPE zC = SCALAR_LOAD((_point_outer + p_inner + 2 * npts));

         SCALAR_TYPE X_PC = SCALAR_SUB(xP, xC);
         SCALAR_TYPE Y_PC = SCALAR_SUB(yP, yC);
         SCALAR_TYPE Z_PC = SCALAR_SUB(zP, zC);

         X_PC = SCALAR_MUL(X_PC, X_PC);
         X_PC = SCALAR_FMA(Y_PC, Y_PC, X_PC);
         X_PC = SCALAR_FMA(Z_PC, Z_PC, X_PC);
         SCALAR_TYPE TVAL = SCALAR_MUL(RHO, X_PC);

         SCALAR_TYPE t00, t01, t02, t03, t04, TVAL_inv_e;

         // Evaluate Boys function
         boys_element<4>(&TVAL, &TVAL_inv_e, &t04, boys_table);

         // Evaluate VRR Buffer
         SCALAR_TYPE t10, t11, t12, t13, t20, t21, t22, t30, t31, t40, tx, ty;

         t03 = SCALAR_MUL(SCALAR_ADD(SCALAR_MUL(TVAL, t04), TVAL_inv_e), SCALAR_SET1(0.28571428571428569843));
         t02 = SCALAR_MUL(SCALAR_ADD(SCALAR_MUL(TVAL, t03), TVAL_inv_e), SCALAR_SET1(0.40000000000000002220));
         t01 = SCALAR_MUL(SCALAR_ADD(SCALAR_MUL(TVAL, t02), TVAL_inv_e), SCALAR_SET1(0.66666666666666662966));
         t00 = SCALAR_MUL(SCALAR_ADD(SCALAR_MUL(TVAL, t01), TVAL_inv_e), SCALAR_SET1(2.00000000000000000000));

         t00 = SCALAR_MUL(eval, t00);
         t01 = SCALAR_MUL(eval, t01);
         t02 = SCALAR_MUL(eval, t02);
         t03 = SCALAR_MUL(eval, t03);
         t04 = SCALAR_MUL(eval, t04);
         t10 = SCALAR_MUL(X_PA, t00);
         t10 = SCALAR_FNMA(X_PC, t01, t10);
         t11 = SCALAR_MUL(X_PA, t01);
         t11 = SCALAR_FNMA(X_PC, t02, t11);
         t12 = SCALAR_MUL(X_PA, t02);
         t12 = SCALAR_FNMA(X_PC, t03, t12);
         t13 = SCALAR_MUL(X_PA, t03);
         t13 = SCALAR_FNMA(X_PC, t04, t13);
         t20 = SCALAR_MUL(X_PA, t10);
         t20 = SCALAR_FNMA(X_PC, t11, t20);
         tx = SCALAR_SUB(t00, t01);
         ty = SCALAR_SET1(0.5 * 1);
         ty = SCALAR_MUL(ty, RHO_INV);
         t20 = SCALAR_FMA(tx, ty, t20);
         t21 = SCALAR_MUL(X_PA, t11);
         t21 = SCALAR_FNMA(X_PC, t12, t21);
         tx = SCALAR_SUB(t01, t02);
         ty = SCALAR_SET1(0.5 * 1);
         ty = SCALAR_MUL(ty, RHO_INV);
         t21 = SCALAR_FMA(tx, ty, t21);
         t22 = SCALAR_MUL(X_PA, t12);
         t22 = SCALAR_FNMA(X_PC, t13, t22);
         tx = SCALAR_SUB(t02, t03);
         ty = SCALAR_SET1(0.5 * 1);
         ty = SCALAR_MUL(ty, RHO_INV);
         t22 = SCALAR_FMA(tx, ty, t22);
         tx = SCALAR_LOAD((temp + 0 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t20);
         SCALAR_STORE((temp + 0 * blockDim.x + threadIdx.x), tx);
         t30 = SCALAR_MUL(X_PA, t20);
         t30 = SCALAR_FNMA(X_PC, t21, t30);
         tx = SCALAR_SUB(t10, t11);
         ty = SCALAR_SET1(0.5 * 2);
         ty = SCALAR_MUL(ty, RHO_INV);
         t30 = SCALAR_FMA(tx, ty, t30);
         t31 = SCALAR_MUL(X_PA, t21);
         t31 = SCALAR_FNMA(X_PC, t22, t31);
         tx = SCALAR_SUB(t11, t12);
         ty = SCALAR_SET1(0.5 * 2);
         ty = SCALAR_MUL(ty, RHO_INV);
         t31 = SCALAR_FMA(tx, ty, t31);
         tx = SCALAR_LOAD((temp + 6 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t30);
         SCALAR_STORE((temp + 6 * blockDim.x + threadIdx.x), tx);
         t40 = SCALAR_MUL(X_PA, t30);
         t40 = SCALAR_FNMA(X_PC, t31, t40);
         tx = SCALAR_SUB(t20, t21);
         ty = SCALAR_SET1(0.5 * 3);
         ty = SCALAR_MUL(ty, RHO_INV);
         t40 = SCALAR_FMA(tx, ty, t40);
         tx = SCALAR_LOAD((temp + 16 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t40);
         SCALAR_STORE((temp + 16 * blockDim.x + threadIdx.x), tx);
         t40 = SCALAR_MUL(Y_PA, t30);
         t40 = SCALAR_FNMA(Y_PC, t31, t40);
         tx = SCALAR_LOAD((temp + 17 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t40);
         SCALAR_STORE((temp + 17 * blockDim.x + threadIdx.x), tx);
         t40 = SCALAR_MUL(Z_PA, t30);
         t40 = SCALAR_FNMA(Z_PC, t31, t40);
         tx = SCALAR_LOAD((temp + 18 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t40);
         SCALAR_STORE((temp + 18 * blockDim.x + threadIdx.x), tx);
         t30 = SCALAR_MUL(Y_PA, t20);
         t30 = SCALAR_FNMA(Y_PC, t21, t30);
         t31 = SCALAR_MUL(Y_PA, t21);
         t31 = SCALAR_FNMA(Y_PC, t22, t31);
         tx = SCALAR_LOAD((temp + 7 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t30);
         SCALAR_STORE((temp + 7 * blockDim.x + threadIdx.x), tx);
         t40 = SCALAR_MUL(Y_PA, t30);
         t40 = SCALAR_FNMA(Y_PC, t31, t40);
         tx = SCALAR_SUB(t20, t21);
         ty = SCALAR_SET1(0.5 * 1);
         ty = SCALAR_MUL(ty, RHO_INV);
         t40 = SCALAR_FMA(tx, ty, t40);
         tx = SCALAR_LOAD((temp + 19 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t40);
         SCALAR_STORE((temp + 19 * blockDim.x + threadIdx.x), tx);
         t40 = SCALAR_MUL(Z_PA, t30);
         t40 = SCALAR_FNMA(Z_PC, t31, t40);
         tx = SCALAR_LOAD((temp + 20 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t40);
         SCALAR_STORE((temp + 20 * blockDim.x + threadIdx.x), tx);
         t30 = SCALAR_MUL(Z_PA, t20);
         t30 = SCALAR_FNMA(Z_PC, t21, t30);
         t31 = SCALAR_MUL(Z_PA, t21);
         t31 = SCALAR_FNMA(Z_PC, t22, t31);
         tx = SCALAR_LOAD((temp + 8 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t30);
         SCALAR_STORE((temp + 8 * blockDim.x + threadIdx.x), tx);
         t40 = SCALAR_MUL(Z_PA, t30);
         t40 = SCALAR_FNMA(Z_PC, t31, t40);
         tx = SCALAR_SUB(t20, t21);
         ty = SCALAR_SET1(0.5 * 1);
         ty = SCALAR_MUL(ty, RHO_INV);
         t40 = SCALAR_FMA(tx, ty, t40);
         tx = SCALAR_LOAD((temp + 21 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t40);
         SCALAR_STORE((temp + 21 * blockDim.x + threadIdx.x), tx);
         t20 = SCALAR_MUL(Y_PA, t10);
         t20 = SCALAR_FNMA(Y_PC, t11, t20);
         t21 = SCALAR_MUL(Y_PA, t11);
         t21 = SCALAR_FNMA(Y_PC, t12, t21);
         t22 = SCALAR_MUL(Y_PA, t12);
         t22 = SCALAR_FNMA(Y_PC, t13, t22);
         tx = SCALAR_LOAD((temp + 1 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t20);
         SCALAR_STORE((temp + 1 * blockDim.x + threadIdx.x), tx);
         t30 = SCALAR_MUL(Y_PA, t20);
         t30 = SCALAR_FNMA(Y_PC, t21, t30);
         tx = SCALAR_SUB(t10, t11);
         ty = SCALAR_SET1(0.5 * 1);
         ty = SCALAR_MUL(ty, RHO_INV);
         t30 = SCALAR_FMA(tx, ty, t30);
         t31 = SCALAR_MUL(Y_PA, t21);
         t31 = SCALAR_FNMA(Y_PC, t22, t31);
         tx = SCALAR_SUB(t11, t12);
         ty = SCALAR_SET1(0.5 * 1);
         ty = SCALAR_MUL(ty, RHO_INV);
         t31 = SCALAR_FMA(tx, ty, t31);
         tx = SCALAR_LOAD((temp + 9 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t30);
         SCALAR_STORE((temp + 9 * blockDim.x + threadIdx.x), tx);
         t40 = SCALAR_MUL(Y_PA, t30);
         t40 = SCALAR_FNMA(Y_PC, t31, t40);
         tx = SCALAR_SUB(t20, t21);
         ty = SCALAR_SET1(0.5 * 2);
         ty = SCALAR_MUL(ty, RHO_INV);
         t40 = SCALAR_FMA(tx, ty, t40);
         tx = SCALAR_LOAD((temp + 22 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t40);
         SCALAR_STORE((temp + 22 * blockDim.x + threadIdx.x), tx);
         t40 = SCALAR_MUL(Z_PA, t30);
         t40 = SCALAR_FNMA(Z_PC, t31, t40);
         tx = SCALAR_LOAD((temp + 23 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t40);
         SCALAR_STORE((temp + 23 * blockDim.x + threadIdx.x), tx);
         t30 = SCALAR_MUL(Z_PA, t20);
         t30 = SCALAR_FNMA(Z_PC, t21, t30);
         t31 = SCALAR_MUL(Z_PA, t21);
         t31 = SCALAR_FNMA(Z_PC, t22, t31);
         tx = SCALAR_LOAD((temp + 10 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t30);
         SCALAR_STORE((temp + 10 * blockDim.x + threadIdx.x), tx);
         t40 = SCALAR_MUL(Z_PA, t30);
         t40 = SCALAR_FNMA(Z_PC, t31, t40);
         tx = SCALAR_SUB(t20, t21);
         ty = SCALAR_SET1(0.5 * 1);
         ty = SCALAR_MUL(ty, RHO_INV);
         t40 = SCALAR_FMA(tx, ty, t40);
         tx = SCALAR_LOAD((temp + 24 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t40);
         SCALAR_STORE((temp + 24 * blockDim.x + threadIdx.x), tx);
         t20 = SCALAR_MUL(Z_PA, t10);
         t20 = SCALAR_FNMA(Z_PC, t11, t20);
         t21 = SCALAR_MUL(Z_PA, t11);
         t21 = SCALAR_FNMA(Z_PC, t12, t21);
         t22 = SCALAR_MUL(Z_PA, t12);
         t22 = SCALAR_FNMA(Z_PC, t13, t22);
         tx = SCALAR_LOAD((temp + 2 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t20);
         SCALAR_STORE((temp + 2 * blockDim.x + threadIdx.x), tx);
         t30 = SCALAR_MUL(Z_PA, t20);
         t30 = SCALAR_FNMA(Z_PC, t21, t30);
         tx = SCALAR_SUB(t10, t11);
         ty = SCALAR_SET1(0.5 * 1);
         ty = SCALAR_MUL(ty, RHO_INV);
         t30 = SCALAR_FMA(tx, ty, t30);
         t31 = SCALAR_MUL(Z_PA, t21);
         t31 = SCALAR_FNMA(Z_PC, t22, t31);
         tx = SCALAR_SUB(t11, t12);
         ty = SCALAR_SET1(0.5 * 1);
         ty = SCALAR_MUL(ty, RHO_INV);
         t31 = SCALAR_FMA(tx, ty, t31);
         tx = SCALAR_LOAD((temp + 11 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t30);
         SCALAR_STORE((temp + 11 * blockDim.x + threadIdx.x), tx);
         t40 = SCALAR_MUL(Z_PA, t30);
         t40 = SCALAR_FNMA(Z_PC, t31, t40);
         tx = SCALAR_SUB(t20, t21);
         ty = SCALAR_SET1(0.5 * 2);
         ty = SCALAR_MUL(ty, RHO_INV);
         t40 = SCALAR_FMA(tx, ty, t40);
         tx = SCALAR_LOAD((temp + 25 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t40);
         SCALAR_STORE((temp + 25 * blockDim.x + threadIdx.x), tx);
         t10 = SCALAR_MUL(Y_PA, t00);
         t10 = SCALAR_FNMA(Y_PC, t01, t10);
         t11 = SCALAR_MUL(Y_PA, t01);
         t11 = SCALAR_FNMA(Y_PC, t02, t11);
         t12 = SCALAR_MUL(Y_PA, t02);
         t12 = SCALAR_FNMA(Y_PC, t03, t12);
         t13 = SCALAR_MUL(Y_PA, t03);
         t13 = SCALAR_FNMA(Y_PC, t04, t13);
         t20 = SCALAR_MUL(Y_PA, t10);
         t20 = SCALAR_FNMA(Y_PC, t11, t20);
         tx = SCALAR_SUB(t00, t01);
         ty = SCALAR_SET1(0.5 * 1);
         ty = SCALAR_MUL(ty, RHO_INV);
         t20 = SCALAR_FMA(tx, ty, t20);
         t21 = SCALAR_MUL(Y_PA, t11);
         t21 = SCALAR_FNMA(Y_PC, t12, t21);
         tx = SCALAR_SUB(t01, t02);
         ty = SCALAR_SET1(0.5 * 1);
         ty = SCALAR_MUL(ty, RHO_INV);
         t21 = SCALAR_FMA(tx, ty, t21);
         t22 = SCALAR_MUL(Y_PA, t12);
         t22 = SCALAR_FNMA(Y_PC, t13, t22);
         tx = SCALAR_SUB(t02, t03);
         ty = SCALAR_SET1(0.5 * 1);
         ty = SCALAR_MUL(ty, RHO_INV);
         t22 = SCALAR_FMA(tx, ty, t22);
         tx = SCALAR_LOAD((temp + 3 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t20);
         SCALAR_STORE((temp + 3 * blockDim.x + threadIdx.x), tx);
         t30 = SCALAR_MUL(Y_PA, t20);
         t30 = SCALAR_FNMA(Y_PC, t21, t30);
         tx = SCALAR_SUB(t10, t11);
         ty = SCALAR_SET1(0.5 * 2);
         ty = SCALAR_MUL(ty, RHO_INV);
         t30 = SCALAR_FMA(tx, ty, t30);
         t31 = SCALAR_MUL(Y_PA, t21);
         t31 = SCALAR_FNMA(Y_PC, t22, t31);
         tx = SCALAR_SUB(t11, t12);
         ty = SCALAR_SET1(0.5 * 2);
         ty = SCALAR_MUL(ty, RHO_INV);
         t31 = SCALAR_FMA(tx, ty, t31);
         tx = SCALAR_LOAD((temp + 12 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t30);
         SCALAR_STORE((temp + 12 * blockDim.x + threadIdx.x), tx);
         t40 = SCALAR_MUL(Y_PA, t30);
         t40 = SCALAR_FNMA(Y_PC, t31, t40);
         tx = SCALAR_SUB(t20, t21);
         ty = SCALAR_SET1(0.5 * 3);
         ty = SCALAR_MUL(ty, RHO_INV);
         t40 = SCALAR_FMA(tx, ty, t40);
         tx = SCALAR_LOAD((temp + 26 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t40);
         SCALAR_STORE((temp + 26 * blockDim.x + threadIdx.x), tx);
         t40 = SCALAR_MUL(Z_PA, t30);
         t40 = SCALAR_FNMA(Z_PC, t31, t40);
         tx = SCALAR_LOAD((temp + 27 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t40);
         SCALAR_STORE((temp + 27 * blockDim.x + threadIdx.x), tx);
         t30 = SCALAR_MUL(Z_PA, t20);
         t30 = SCALAR_FNMA(Z_PC, t21, t30);
         t31 = SCALAR_MUL(Z_PA, t21);
         t31 = SCALAR_FNMA(Z_PC, t22, t31);
         tx = SCALAR_LOAD((temp + 13 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t30);
         SCALAR_STORE((temp + 13 * blockDim.x + threadIdx.x), tx);
         t40 = SCALAR_MUL(Z_PA, t30);
         t40 = SCALAR_FNMA(Z_PC, t31, t40);
         tx = SCALAR_SUB(t20, t21);
         ty = SCALAR_SET1(0.5 * 1);
         ty = SCALAR_MUL(ty, RHO_INV);
         t40 = SCALAR_FMA(tx, ty, t40);
         tx = SCALAR_LOAD((temp + 28 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t40);
         SCALAR_STORE((temp + 28 * blockDim.x + threadIdx.x), tx);
         t20 = SCALAR_MUL(Z_PA, t10);
         t20 = SCALAR_FNMA(Z_PC, t11, t20);
         t21 = SCALAR_MUL(Z_PA, t11);
         t21 = SCALAR_FNMA(Z_PC, t12, t21);
         t22 = SCALAR_MUL(Z_PA, t12);
         t22 = SCALAR_FNMA(Z_PC, t13, t22);
         tx = SCALAR_LOAD((temp + 4 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t20);
         SCALAR_STORE((temp + 4 * blockDim.x + threadIdx.x), tx);
         t30 = SCALAR_MUL(Z_PA, t20);
         t30 = SCALAR_FNMA(Z_PC, t21, t30);
         tx = SCALAR_SUB(t10, t11);
         ty = SCALAR_SET1(0.5 * 1);
         ty = SCALAR_MUL(ty, RHO_INV);
         t30 = SCALAR_FMA(tx, ty, t30);
         t31 = SCALAR_MUL(Z_PA, t21);
         t31 = SCALAR_FNMA(Z_PC, t22, t31);
         tx = SCALAR_SUB(t11, t12);
         ty = SCALAR_SET1(0.5 * 1);
         ty = SCALAR_MUL(ty, RHO_INV);
         t31 = SCALAR_FMA(tx, ty, t31);
         tx = SCALAR_LOAD((temp + 14 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t30);
         SCALAR_STORE((temp + 14 * blockDim.x + threadIdx.x), tx);
         t40 = SCALAR_MUL(Z_PA, t30);
         t40 = SCALAR_FNMA(Z_PC, t31, t40);
         tx = SCALAR_SUB(t20, t21);
         ty = SCALAR_SET1(0.5 * 2);
         ty = SCALAR_MUL(ty, RHO_INV);
         t40 = SCALAR_FMA(tx, ty, t40);
         tx = SCALAR_LOAD((temp + 29 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t40);
         SCALAR_STORE((temp + 29 * blockDim.x + threadIdx.x), tx);
         t10 = SCALAR_MUL(Z_PA, t00);
         t10 = SCALAR_FNMA(Z_PC, t01, t10);
         t11 = SCALAR_MUL(Z_PA, t01);
         t11 = SCALAR_FNMA(Z_PC, t02, t11);
         t12 = SCALAR_MUL(Z_PA, t02);
         t12 = SCALAR_FNMA(Z_PC, t03, t12);
         t13 = SCALAR_MUL(Z_PA, t03);
         t13 = SCALAR_FNMA(Z_PC, t04, t13);
         t20 = SCALAR_MUL(Z_PA, t10);
         t20 = SCALAR_FNMA(Z_PC, t11, t20);
         tx = SCALAR_SUB(t00, t01);
         ty = SCALAR_SET1(0.5 * 1);
         ty = SCALAR_MUL(ty, RHO_INV);
         t20 = SCALAR_FMA(tx, ty, t20);
         t21 = SCALAR_MUL(Z_PA, t11);
         t21 = SCALAR_FNMA(Z_PC, t12, t21);
         tx = SCALAR_SUB(t01, t02);
         ty = SCALAR_SET1(0.5 * 1);
         ty = SCALAR_MUL(ty, RHO_INV);
         t21 = SCALAR_FMA(tx, ty, t21);
         t22 = SCALAR_MUL(Z_PA, t12);
         t22 = SCALAR_FNMA(Z_PC, t13, t22);
         tx = SCALAR_SUB(t02, t03);
         ty = SCALAR_SET1(0.5 * 1);
         ty = SCALAR_MUL(ty, RHO_INV);
         t22 = SCALAR_FMA(tx, ty, t22);
         tx = SCALAR_LOAD((temp + 5 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t20);
         SCALAR_STORE((temp + 5 * blockDim.x + threadIdx.x), tx);
         t30 = SCALAR_MUL(Z_PA, t20);
         t30 = SCALAR_FNMA(Z_PC, t21, t30);
         tx = SCALAR_SUB(t10, t11);
         ty = SCALAR_SET1(0.5 * 2);
         ty = SCALAR_MUL(ty, RHO_INV);
         t30 = SCALAR_FMA(tx, ty, t30);
         t31 = SCALAR_MUL(Z_PA, t21);
         t31 = SCALAR_FNMA(Z_PC, t22, t31);
         tx = SCALAR_SUB(t11, t12);
         ty = SCALAR_SET1(0.5 * 2);
         ty = SCALAR_MUL(ty, RHO_INV);
         t31 = SCALAR_FMA(tx, ty, t31);
         tx = SCALAR_LOAD((temp + 15 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t30);
         SCALAR_STORE((temp + 15 * blockDim.x + threadIdx.x), tx);
         t40 = SCALAR_MUL(Z_PA, t30);
         t40 = SCALAR_FNMA(Z_PC, t31, t40);
         tx = SCALAR_SUB(t20, t21);
         ty = SCALAR_SET1(0.5 * 3);
         ty = SCALAR_MUL(ty, RHO_INV);
         t40 = SCALAR_FMA(tx, ty, t40);
         tx = SCALAR_LOAD((temp + 30 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t40);
         SCALAR_STORE((temp + 30 * blockDim.x + threadIdx.x), tx);
      }

      double *Xik = (Xi + p_outer + p_inner);
      double *Xjk = (Xj + p_outer + p_inner);
      double *Gik = (Gi + p_outer + p_inner);
      double *Gjk = (Gj + p_outer + p_inner);

      SCALAR_TYPE const_value_v = SCALAR_LOAD((weights + p_outer + p_inner));

      for(int c0 = 0; c0 <= 2; ++c0) {
         for(int c1 = 0; c1 <= c0; ++c1) {
            int m = 2 - c0;
            int n = c0 - c1;
            int p = c1;

            int idxB = (((2 - m) * (2 - m + 1)) >> 1) + p;

            double X_ABp = 1.0, comb_m_i = 1.0;
            for(int i = 0; i <= m; ++i) {
               double rcp_i;

               double Y_ABp = 1.0, comb_n_j = 1.0;
               for(int j = 0; j <= n; ++j) {
                  double rcp_j;

                  double Z_ABp = 1.0, comb_p_k = 1.0;
                  for(int k = 0; k <= p; ++k) {
                     double rcp_k;
                     int mv, pv, Lv = 4 - i - j - k;

                     int offset = (Lv * (Lv + 1) * (Lv + 2) - 24) / 6;
                     double const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
                     SCALAR_TYPE tx, ty, tz, tw;
                     SCALAR_TYPE const_value_w = SCALAR_MUL(const_value_v, const_value);

                     mv = 2 + m - i; pv = 0 + p - k;
                     tx = SCALAR_LOAD((Xik + 0 * ldX));
                     ty = SCALAR_LOAD((Xjk + idxB * ldX));
                     tz = SCALAR_LOAD((Gik + 0 * ldG));
                     tw = SCALAR_LOAD((Gjk + idxB * ldG));
                     SCALAR_TYPE t0 = SCALAR_LOAD((temp + (offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * blockDim.x + threadIdx.x));
                     t0 = SCALAR_MUL(t0, const_value_w);
                     tz = SCALAR_FMA(ty, t0, tz);
                     tw = SCALAR_FMA(tx, t0, tw);
                     SCALAR_STORE((Gik + 0 * ldG), tz);
                     SCALAR_STORE((Gjk + idxB * ldG), tw);
                     mv = 1 + m - i; pv = 0 + p - k;
                     tx = SCALAR_LOAD((Xik + 1 * ldX));
                     ty = SCALAR_LOAD((Xjk + idxB * ldX));
                     tz = SCALAR_LOAD((Gik + 1 * ldG));
                     tw = SCALAR_LOAD((Gjk + idxB * ldG));
                     SCALAR_TYPE t1 = SCALAR_LOAD((temp + (offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * blockDim.x + threadIdx.x));
                     t1 = SCALAR_MUL(t1, const_value_w);
                     tz = SCALAR_FMA(ty, t1, tz);
                     tw = SCALAR_FMA(tx, t1, tw);
                     SCALAR_STORE((Gik + 1 * ldG), tz);
                     SCALAR_STORE((Gjk + idxB * ldG), tw);
                     mv = 1 + m - i; pv = 1 + p - k;
                     tx = SCALAR_LOAD((Xik + 2 * ldX));
                     ty = SCALAR_LOAD((Xjk + idxB * ldX));
                     tz = SCALAR_LOAD((Gik + 2 * ldG));
                     tw = SCALAR_LOAD((Gjk + idxB * ldG));
                     SCALAR_TYPE t2 = SCALAR_LOAD((temp + (offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * blockDim.x + threadIdx.x));
                     t2 = SCALAR_MUL(t2, const_value_w);
                     tz = SCALAR_FMA(ty, t2, tz);
                     tw = SCALAR_FMA(tx, t2, tw);
                     SCALAR_STORE((Gik + 2 * ldG), tz);
                     SCALAR_STORE((Gjk + idxB * ldG), tw);
                     mv = 0 + m - i; pv = 0 + p - k;
                     tx = SCALAR_LOAD((Xik + 3 * ldX));
                     ty = SCALAR_LOAD((Xjk + idxB * ldX));
                     tz = SCALAR_LOAD((Gik + 3 * ldG));
                     tw = SCALAR_LOAD((Gjk + idxB * ldG));
                     SCALAR_TYPE t3 = SCALAR_LOAD((temp + (offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * blockDim.x + threadIdx.x));
                     t3 = SCALAR_MUL(t3, const_value_w);
                     tz = SCALAR_FMA(ty, t3, tz);
                     tw = SCALAR_FMA(tx, t3, tw);
                     SCALAR_STORE((Gik + 3 * ldG), tz);
                     SCALAR_STORE((Gjk + idxB * ldG), tw);
                     mv = 0 + m - i; pv = 1 + p - k;
                     tx = SCALAR_LOAD((Xik + 4 * ldX));
                     ty = SCALAR_LOAD((Xjk + idxB * ldX));
                     tz = SCALAR_LOAD((Gik + 4 * ldG));
                     tw = SCALAR_LOAD((Gjk + idxB * ldG));
                     SCALAR_TYPE t4 = SCALAR_LOAD((temp + (offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * blockDim.x + threadIdx.x));
                     t4 = SCALAR_MUL(t4, const_value_w);
                     tz = SCALAR_FMA(ty, t4, tz);
                     tw = SCALAR_FMA(tx, t4, tw);
                     SCALAR_STORE((Gik + 4 * ldG), tz);
                     SCALAR_STORE((Gjk + idxB * ldG), tw);
                     mv = 0 + m - i; pv = 2 + p - k;
                     tx = SCALAR_LOAD((Xik + 5 * ldX));
                     ty = SCALAR_LOAD((Xjk + idxB * ldX));
                     tz = SCALAR_LOAD((Gik + 5 * ldG));
                     tw = SCALAR_LOAD((Gjk + idxB * ldG));
                     SCALAR_TYPE t5 = SCALAR_LOAD((temp + (offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * blockDim.x + threadIdx.x));
                     t5 = SCALAR_MUL(t5, const_value_w);
                     tz = SCALAR_FMA(ty, t5, tz);
                     tw = SCALAR_FMA(tx, t5, tw);
                     SCALAR_STORE((Gik + 5 * ldG), tz);
                     SCALAR_STORE((Gjk + idxB * ldG), tw);

                     Z_ABp = SCALAR_MUL(Z_ABp, Z_AB);
                     rcp_k = SCALAR_RECIPROCAL(k + 1);
                     comb_p_k = SCALAR_MUL(comb_p_k, p - k);
                     comb_p_k = SCALAR_MUL(comb_p_k, rcp_k);
                  }

                  Y_ABp = SCALAR_MUL(Y_ABp, Y_AB);
                  rcp_j = SCALAR_RECIPROCAL(j + 1);
                  comb_n_j = SCALAR_MUL(comb_n_j, n - j);
                  comb_n_j = SCALAR_MUL(comb_n_j, rcp_j);
               }

               X_ABp = SCALAR_MUL(X_ABp, X_AB);
               rcp_i = SCALAR_RECIPROCAL(i + 1);
               comb_m_i = SCALAR_MUL(comb_m_i, m - i);
               comb_m_i = SCALAR_MUL(comb_m_i, rcp_i);
            }
         }
      }
   }
}
}
