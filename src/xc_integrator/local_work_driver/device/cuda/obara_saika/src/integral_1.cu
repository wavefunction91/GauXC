#include <math.h>
#include "chebyshev_boys_computation.hpp"
#include "integral_data_types.hpp"
#include "config_obara_saika.hpp"
#include "integral_1.hu"

#define PI 3.14159265358979323846

#define MIN(a,b)			\
  ({ __typeof__ (a) _a = (a);	        \
  __typeof__ (b) _b = (b);		\
  _a < _b ? _a : _b; })

namespace XGPU {
__global__ void integral_1(size_t npts,
                          point rA,
                          point rB,
                          point rAB,
                          int nprim_pair,
                          prim_pair *ppair,
                          double *_points,
                          double *Xi,
                          int ldX,
                          double *Gi,
                          int ldG, 
                          double *weights,
                          double *boys_table) {
   __shared__ double *temp;
   for(size_t p_outer = blockIdx.x * blockDim.x; p_outer < npts; p_outer += gridDim.x * blockDim.x) {
      double *_point_outer = (_points + p_outer);

      size_t p_inner = (threadIdx.x < (npts - p_outer)) ? threadIdx.x : (npts - p_outer);

      double xA = shpair.rA.x;
      double yA = shpair.rA.y;
      double zA = shpair.rA.z;

      for(int i = 0; i < 9; ++i) SCALAR_STORE((temp + i * blockDim.x + threadIdx.x), SCALAR_ZERO());

      for(int ij = 0; ij < shpair.nprim_pair; ++ij) {
         double RHO = shpair.prim_pairs[ij].gamma;
         double RHO_INV = 1.0 / RHO;

         constexpr double X_PA = 0.0;
         constexpr double Y_PA = 0.0;
         constexpr double Z_PA = 0.0;

         double eval = shpair.prim_pairs[ij].K_coeff_prod;

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

         SCALAR_TYPE t00, t01, t02, TVAL_inv_e;

         // Evaluate Boys function
         boys_element<2>(&TVAL, &TVAL_inv_e, &t02, boys_table);

         // Evaluate VRR Buffer
         SCALAR_TYPE t10, t11, t20, tx, ty;

         t01 = SCALAR_MUL(SCALAR_ADD(SCALAR_MUL(TVAL, t02), TVAL_inv_e), SCALAR_SET1(0.66666666666666662966));
         t00 = SCALAR_MUL(SCALAR_ADD(SCALAR_MUL(TVAL, t01), TVAL_inv_e), SCALAR_SET1(2.00000000000000000000));

         t00 = SCALAR_MUL(eval, t00);
         t01 = SCALAR_MUL(eval, t01);
         t02 = SCALAR_MUL(eval, t02);
         t10 = SCALAR_MUL(X_PA, t00);
         t10 = SCALAR_FNMA(X_PC, t01, t10);
         t11 = SCALAR_MUL(X_PA, t01);
         t11 = SCALAR_FNMA(X_PC, t02, t11);
         tx = SCALAR_LOAD((temp + 0 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t10);
         SCALAR_STORE((temp + 0 * blockDim.x + threadIdx.x), tx);
         t20 = SCALAR_MUL(X_PA, t10);
         t20 = SCALAR_FNMA(X_PC, t11, t20);
         tx = SCALAR_SUB(t00, t01);
         ty = SCALAR_SET1(0.5 * 1);
         ty = SCALAR_MUL(ty, RHO_INV);
         t20 = SCALAR_FMA(tx, ty, t20);
         tx = SCALAR_LOAD((temp + 3 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t20);
         SCALAR_STORE((temp + 3 * blockDim.x + threadIdx.x), tx);
         t20 = SCALAR_MUL(Y_PA, t10);
         t20 = SCALAR_FNMA(Y_PC, t11, t20);
         tx = SCALAR_LOAD((temp + 4 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t20);
         SCALAR_STORE((temp + 4 * blockDim.x + threadIdx.x), tx);
         t20 = SCALAR_MUL(Z_PA, t10);
         t20 = SCALAR_FNMA(Z_PC, t11, t20);
         tx = SCALAR_LOAD((temp + 5 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t20);
         SCALAR_STORE((temp + 5 * blockDim.x + threadIdx.x), tx);
         t10 = SCALAR_MUL(Y_PA, t00);
         t10 = SCALAR_FNMA(Y_PC, t01, t10);
         t11 = SCALAR_MUL(Y_PA, t01);
         t11 = SCALAR_FNMA(Y_PC, t02, t11);
         tx = SCALAR_LOAD((temp + 1 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t10);
         SCALAR_STORE((temp + 1 * blockDim.x + threadIdx.x), tx);
         t20 = SCALAR_MUL(Y_PA, t10);
         t20 = SCALAR_FNMA(Y_PC, t11, t20);
         tx = SCALAR_SUB(t00, t01);
         ty = SCALAR_SET1(0.5 * 1);
         ty = SCALAR_MUL(ty, RHO_INV);
         t20 = SCALAR_FMA(tx, ty, t20);
         tx = SCALAR_LOAD((temp + 6 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t20);
         SCALAR_STORE((temp + 6 * blockDim.x + threadIdx.x), tx);
         t20 = SCALAR_MUL(Z_PA, t10);
         t20 = SCALAR_FNMA(Z_PC, t11, t20);
         tx = SCALAR_LOAD((temp + 7 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t20);
         SCALAR_STORE((temp + 7 * blockDim.x + threadIdx.x), tx);
         t10 = SCALAR_MUL(Z_PA, t00);
         t10 = SCALAR_FNMA(Z_PC, t01, t10);
         t11 = SCALAR_MUL(Z_PA, t01);
         t11 = SCALAR_FNMA(Z_PC, t02, t11);
         tx = SCALAR_LOAD((temp + 2 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t10);
         SCALAR_STORE((temp + 2 * blockDim.x + threadIdx.x), tx);
         t20 = SCALAR_MUL(Z_PA, t10);
         t20 = SCALAR_FNMA(Z_PC, t11, t20);
         tx = SCALAR_SUB(t00, t01);
         ty = SCALAR_SET1(0.5 * 1);
         ty = SCALAR_MUL(ty, RHO_INV);
         t20 = SCALAR_FMA(tx, ty, t20);
         tx = SCALAR_LOAD((temp + 8 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t20);
         SCALAR_STORE((temp + 8 * blockDim.x + threadIdx.x), tx);
      }

      double *Xik = (Xi + p_outer + p_inner);
      double *Gik = (Gi + p_outer + p_inner);

      SCALAR_TYPE tx, wg, xik, gik;
      tx  = SCALAR_LOAD((temp + 3 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 0 * ldX));
      gik = SCALAR_LOAD((Gik + 0 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 0 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 4 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 0 * ldX));
      gik = SCALAR_LOAD((Gik + 1 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 1 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 5 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 0 * ldX));
      gik = SCALAR_LOAD((Gik + 2 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 2 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 4 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 1 * ldX));
      gik = SCALAR_LOAD((Gik + 0 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 0 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 6 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 1 * ldX));
      gik = SCALAR_LOAD((Gik + 1 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 1 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 7 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 1 * ldX));
      gik = SCALAR_LOAD((Gik + 2 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 2 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 5 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 2 * ldX));
      gik = SCALAR_LOAD((Gik + 0 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 0 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 7 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 2 * ldX));
      gik = SCALAR_LOAD((Gik + 1 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 1 * ldG), gik);
      tx  = SCALAR_LOAD((temp + 8 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 2 * ldX));
      gik = SCALAR_LOAD((Gik + 2 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 2 * ldG), gik);
   }
}
}
