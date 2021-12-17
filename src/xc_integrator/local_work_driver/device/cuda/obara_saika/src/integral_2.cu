#include <math.h>
#include "chebyshev_boys_computation.hpp"
#include "integral_data_types.hpp"
#include "config_obara_saika.hpp"
#include "integral_2.hu"

#define PI 3.14159265358979323846

#define MIN(a,b)			\
  ({ __typeof__ (a) _a = (a);	        \
  __typeof__ (b) _b = (b);		\
  _a < _b ? _a : _b; })

__global__ void integral_2(size_t npts,
                          shell_pair *shpair,
                          double *_points,
                          double *Xi,
                          int ldX,
                          double *Gi,
                          int ldG, 
                          double *weights) {
   __shared__ double *temp;
   for(size_t p_outer = blockIdx.x * blockDim.x; p_outer < npts; p_outer += gridDim.x * blockDim.x) {
      double *_point_outer = (_points + p_outer);

      size_t p_inner = (threadIdx.x < (npts - p_outer)) ? threadIdx.x : (npts - p_outer);

      double xA = shpair[0].rA.x;
      double yA = shpair[0].rA.y;
      double zA = shpair[0].rA.z;

      for(int i = 0; i < 31; ++i) SCALAR_STORE((temp + i * blockDim.x + threadIdx.x), SCALAR_ZERO());

      for(int ij = 0; ij < shpair[0].nprim_pair; ++ij) {
         double RHO = shpair[0].prim_pairs[ij].gamma;
         double RHO_INV = 1.0 / RHO;

         constexpr double X_PA = 0.0;
         constexpr double Y_PA = 0.0;
         constexpr double Z_PA = 0.0;

         double eval = shpair[0].prim_pairs[ij].coeff_prod * shpair[0].prim_pairs[ij].K;

         // Evaluate T Values
         SCALAR_TYPE xC = SCALAR_LOAD((_point_outer + p_inner + 0 * npts));
         SCALAR_TYPE yC = SCALAR_LOAD((_point_outer + p_inner + 1 * npts));
         SCALAR_TYPE zC = SCALAR_LOAD((_point_outer + p_inner + 2 * npts));

         SCALAR_TYPE X_PC = SCALAR_SUB(xA, xC);
         SCALAR_TYPE Y_PC = SCALAR_SUB(yA, yC);
         SCALAR_TYPE Z_PC = SCALAR_SUB(zA, zC);

         X_PC = SCALAR_MUL(X_PC, X_PC);
         X_PC = SCALAR_FMA(Y_PC, Y_PC, X_PC);
         X_PC = SCALAR_FMA(Z_PC, Z_PC, X_PC);
         SCALAR_TYPE TVAL = SCALAR_MUL(RHO, X_PC);

         SCALAR_TYPE t00, t01, t02, t03, t04;

         // Evaluate Boys function
         t00 = GauXC::gauxc_boys_element<0>(TVAL);
         t01 = GauXC::gauxc_boys_element<1>(TVAL);
         t02 = GauXC::gauxc_boys_element<2>(TVAL);
         t03 = GauXC::gauxc_boys_element<3>(TVAL);
         t04 = GauXC::gauxc_boys_element<4>(TVAL);

         // Evaluate VRR Buffer
         SCALAR_TYPE t10, t11, t12, t13, t20, t21, t22, t30, t31, t40, tx, ty;

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
      double *Gik = (Gi + p_outer + p_inner);

      SCALAR_TYPE tx, wg, xik, gik;
      tx  = SCALAR_LOAD((temp + 16 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 0 * ldX));
      gik = SCALAR_LOAD((Gik + 0 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 0 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 17 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 0 * ldX));
      gik = SCALAR_LOAD((Gik + 1 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 1 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 18 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 0 * ldX));
      gik = SCALAR_LOAD((Gik + 2 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 2 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 19 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 0 * ldX));
      gik = SCALAR_LOAD((Gik + 3 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 3 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 20 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 0 * ldX));
      gik = SCALAR_LOAD((Gik + 4 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 4 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 21 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 0 * ldX));
      gik = SCALAR_LOAD((Gik + 5 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 5 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 17 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 1 * ldX));
      gik = SCALAR_LOAD((Gik + 0 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 0 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 19 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 1 * ldX));
      gik = SCALAR_LOAD((Gik + 1 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 1 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 20 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 1 * ldX));
      gik = SCALAR_LOAD((Gik + 2 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 2 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 22 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 1 * ldX));
      gik = SCALAR_LOAD((Gik + 3 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 3 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 23 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 1 * ldX));
      gik = SCALAR_LOAD((Gik + 4 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 4 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 24 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 1 * ldX));
      gik = SCALAR_LOAD((Gik + 5 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 5 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 18 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 2 * ldX));
      gik = SCALAR_LOAD((Gik + 0 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 0 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 20 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 2 * ldX));
      gik = SCALAR_LOAD((Gik + 1 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 1 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 21 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 2 * ldX));
      gik = SCALAR_LOAD((Gik + 2 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 2 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 23 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 2 * ldX));
      gik = SCALAR_LOAD((Gik + 3 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 3 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 24 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 2 * ldX));
      gik = SCALAR_LOAD((Gik + 4 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 4 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 25 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 2 * ldX));
      gik = SCALAR_LOAD((Gik + 5 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 5 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 19 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 3 * ldX));
      gik = SCALAR_LOAD((Gik + 0 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 0 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 22 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 3 * ldX));
      gik = SCALAR_LOAD((Gik + 1 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 1 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 23 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 3 * ldX));
      gik = SCALAR_LOAD((Gik + 2 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 2 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 26 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 3 * ldX));
      gik = SCALAR_LOAD((Gik + 3 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 3 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 27 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 3 * ldX));
      gik = SCALAR_LOAD((Gik + 4 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 4 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 28 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 3 * ldX));
      gik = SCALAR_LOAD((Gik + 5 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 5 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 20 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 4 * ldX));
      gik = SCALAR_LOAD((Gik + 0 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 0 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 23 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 4 * ldX));
      gik = SCALAR_LOAD((Gik + 1 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 1 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 24 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 4 * ldX));
      gik = SCALAR_LOAD((Gik + 2 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 2 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 27 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 4 * ldX));
      gik = SCALAR_LOAD((Gik + 3 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 3 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 28 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 4 * ldX));
      gik = SCALAR_LOAD((Gik + 4 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 4 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 29 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 4 * ldX));
      gik = SCALAR_LOAD((Gik + 5 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 5 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 21 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 5 * ldX));
      gik = SCALAR_LOAD((Gik + 0 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 0 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 24 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 5 * ldX));
      gik = SCALAR_LOAD((Gik + 1 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 1 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 25 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 5 * ldX));
      gik = SCALAR_LOAD((Gik + 2 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 2 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 28 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 5 * ldX));
      gik = SCALAR_LOAD((Gik + 3 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 3 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 29 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 5 * ldX));
      gik = SCALAR_LOAD((Gik + 4 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 4 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 30 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 5 * ldX));
      gik = SCALAR_LOAD((Gik + 5 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 5 * ldG), gik);
   }
}