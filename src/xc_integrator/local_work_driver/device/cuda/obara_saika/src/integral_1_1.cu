#include <math.h>
#include "chebyshev_boys_computation.hpp"
#include "integral_data_types.hpp"
#include "config_obara_saika.hpp"
#include "integral_1_1.hu"

#define PI 3.14159265358979323846

#define MIN(a,b)			\
  ({ __typeof__ (a) _a = (a);	        \
  __typeof__ (b) _b = (b);		\
  _a < _b ? _a : _b; })

__global__ void integral_1_1(size_t npts,
                  shell_pair shpair,
                             double *_points,
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

      double X_AB = shpair.rAB.x;
      double Y_AB = shpair.rAB.y;
      double Z_AB = shpair.rAB.z;

      for(int i = 0; i < 9; ++i) SCALAR_STORE((temp + i * blockDim.x + threadIdx.x), SCALAR_ZERO());

      for(int ij = 0; ij < shpair.nprim_pair; ++ij) {
         double RHO = shpair.prim_pairs[ij].gamma;
         double RHO_INV = 1.0 / RHO;
         double X_PA = shpair.prim_pairs[ij].PA.x;
         double Y_PA = shpair.prim_pairs[ij].PA.y;
         double Z_PA = shpair.prim_pairs[ij].PA.z;

         double xP = shpair.prim_pairs[ij].P.x;
         double yP = shpair.prim_pairs[ij].P.y;
         double zP = shpair.prim_pairs[ij].P.z;

         double eval = shpair.prim_pairs[ij].coeff_prod * shpair.prim_pairs[ij].K;

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

         SCALAR_TYPE t00, t01, t02;

         // Evaluate Boys function
         t00 = GauXC::gauxc_boys_element<0>(boys_table, TVAL);
         t01 = GauXC::gauxc_boys_element<1>(boys_table, TVAL);
         t02 = GauXC::gauxc_boys_element<2>(boys_table, TVAL);

         // Evaluate VRR Buffer
         SCALAR_TYPE t10, t11, t20, tx, ty;

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
      double *Xjk = (Xj + p_outer + p_inner);
      double *Gik = (Gi + p_outer + p_inner);
      double *Gjk = (Gj + p_outer + p_inner);

      SCALAR_TYPE const_value_v = SCALAR_LOAD((weights + p_outer + p_inner));

      double const_value, X_ABp, Y_ABp, Z_ABp, comb_m_i, comb_n_j, comb_p_k;
      SCALAR_TYPE const_value_w;
      SCALAR_TYPE tx, ty, tz, tw, t0, t1, t2;

      X_ABp = 1.0; comb_m_i = 1.0;
      Y_ABp = 1.0; comb_n_j = 1.0;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
      const_value_w = SCALAR_MUL(const_value_v, const_value);
      tx = SCALAR_LOAD((Xik + 0 * ldX));
      ty = SCALAR_LOAD((Xjk + 0 * ldX));
      tz = SCALAR_LOAD((Gik + 0 * ldG));
      tw = SCALAR_LOAD((Gjk + 0 * ldG));
      t0 = SCALAR_LOAD((temp + 3 * blockDim.x + threadIdx.x));
      t0 = SCALAR_MUL(t0, const_value_w);
      tz = SCALAR_FMA(ty, t0, tz);
      tw = SCALAR_FMA(tx, t0, tw);
      SCALAR_STORE((Gik + 0 * ldG), tz);
      SCALAR_STORE((Gjk + 0 * ldG), tw);
      tx = SCALAR_LOAD((Xik + 1 * ldX));
      ty = SCALAR_LOAD((Xjk + 0 * ldX));
      tz = SCALAR_LOAD((Gik + 1 * ldG));
      tw = SCALAR_LOAD((Gjk + 0 * ldG));
      t1 = SCALAR_LOAD((temp + 4 * blockDim.x + threadIdx.x));
      t1 = SCALAR_MUL(t1, const_value_w);
      tz = SCALAR_FMA(ty, t1, tz);
      tw = SCALAR_FMA(tx, t1, tw);
      SCALAR_STORE((Gik + 1 * ldG), tz);
      SCALAR_STORE((Gjk + 0 * ldG), tw);
      tx = SCALAR_LOAD((Xik + 2 * ldX));
      ty = SCALAR_LOAD((Xjk + 0 * ldX));
      tz = SCALAR_LOAD((Gik + 2 * ldG));
      tw = SCALAR_LOAD((Gjk + 0 * ldG));
      t2 = SCALAR_LOAD((temp + 5 * blockDim.x + threadIdx.x));
      t2 = SCALAR_MUL(t2, const_value_w);
      tz = SCALAR_FMA(ty, t2, tz);
      tw = SCALAR_FMA(tx, t2, tw);
      SCALAR_STORE((Gik + 2 * ldG), tz);
      SCALAR_STORE((Gjk + 0 * ldG), tw);
      X_ABp = SCALAR_MUL(X_ABp, X_AB); comb_m_i = SCALAR_MUL(comb_m_i * 1, SCALAR_RECIPROCAL(1));
      Y_ABp = 1.0; comb_n_j = 1.0;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
      const_value_w = SCALAR_MUL(const_value_v, const_value);
      tx = SCALAR_LOAD((Xik + 0 * ldX));
      ty = SCALAR_LOAD((Xjk + 0 * ldX));
      tz = SCALAR_LOAD((Gik + 0 * ldG));
      tw = SCALAR_LOAD((Gjk + 0 * ldG));
      t0 = SCALAR_LOAD((temp + 0 * blockDim.x + threadIdx.x));
      t0 = SCALAR_MUL(t0, const_value_w);
      tz = SCALAR_FMA(ty, t0, tz);
      tw = SCALAR_FMA(tx, t0, tw);
      SCALAR_STORE((Gik + 0 * ldG), tz);
      SCALAR_STORE((Gjk + 0 * ldG), tw);
      tx = SCALAR_LOAD((Xik + 1 * ldX));
      ty = SCALAR_LOAD((Xjk + 0 * ldX));
      tz = SCALAR_LOAD((Gik + 1 * ldG));
      tw = SCALAR_LOAD((Gjk + 0 * ldG));
      t1 = SCALAR_LOAD((temp + 1 * blockDim.x + threadIdx.x));
      t1 = SCALAR_MUL(t1, const_value_w);
      tz = SCALAR_FMA(ty, t1, tz);
      tw = SCALAR_FMA(tx, t1, tw);
      SCALAR_STORE((Gik + 1 * ldG), tz);
      SCALAR_STORE((Gjk + 0 * ldG), tw);
      tx = SCALAR_LOAD((Xik + 2 * ldX));
      ty = SCALAR_LOAD((Xjk + 0 * ldX));
      tz = SCALAR_LOAD((Gik + 2 * ldG));
      tw = SCALAR_LOAD((Gjk + 0 * ldG));
      t2 = SCALAR_LOAD((temp + 2 * blockDim.x + threadIdx.x));
      t2 = SCALAR_MUL(t2, const_value_w);
      tz = SCALAR_FMA(ty, t2, tz);
      tw = SCALAR_FMA(tx, t2, tw);
      SCALAR_STORE((Gik + 2 * ldG), tz);
      SCALAR_STORE((Gjk + 0 * ldG), tw);
      X_ABp = 1.0; comb_m_i = 1.0;
      Y_ABp = 1.0; comb_n_j = 1.0;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
      const_value_w = SCALAR_MUL(const_value_v, const_value);
      tx = SCALAR_LOAD((Xik + 0 * ldX));
      ty = SCALAR_LOAD((Xjk + 1 * ldX));
      tz = SCALAR_LOAD((Gik + 0 * ldG));
      tw = SCALAR_LOAD((Gjk + 1 * ldG));
      t0 = SCALAR_LOAD((temp + 4 * blockDim.x + threadIdx.x));
      t0 = SCALAR_MUL(t0, const_value_w);
      tz = SCALAR_FMA(ty, t0, tz);
      tw = SCALAR_FMA(tx, t0, tw);
      SCALAR_STORE((Gik + 0 * ldG), tz);
      SCALAR_STORE((Gjk + 1 * ldG), tw);
      tx = SCALAR_LOAD((Xik + 1 * ldX));
      ty = SCALAR_LOAD((Xjk + 1 * ldX));
      tz = SCALAR_LOAD((Gik + 1 * ldG));
      tw = SCALAR_LOAD((Gjk + 1 * ldG));
      t1 = SCALAR_LOAD((temp + 6 * blockDim.x + threadIdx.x));
      t1 = SCALAR_MUL(t1, const_value_w);
      tz = SCALAR_FMA(ty, t1, tz);
      tw = SCALAR_FMA(tx, t1, tw);
      SCALAR_STORE((Gik + 1 * ldG), tz);
      SCALAR_STORE((Gjk + 1 * ldG), tw);
      tx = SCALAR_LOAD((Xik + 2 * ldX));
      ty = SCALAR_LOAD((Xjk + 1 * ldX));
      tz = SCALAR_LOAD((Gik + 2 * ldG));
      tw = SCALAR_LOAD((Gjk + 1 * ldG));
      t2 = SCALAR_LOAD((temp + 7 * blockDim.x + threadIdx.x));
      t2 = SCALAR_MUL(t2, const_value_w);
      tz = SCALAR_FMA(ty, t2, tz);
      tw = SCALAR_FMA(tx, t2, tw);
      SCALAR_STORE((Gik + 2 * ldG), tz);
      SCALAR_STORE((Gjk + 1 * ldG), tw);
      Y_ABp = SCALAR_MUL(Y_ABp, Y_AB); comb_n_j = SCALAR_MUL(comb_n_j * 1, SCALAR_RECIPROCAL(1));
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
      const_value_w = SCALAR_MUL(const_value_v, const_value);
      tx = SCALAR_LOAD((Xik + 0 * ldX));
      ty = SCALAR_LOAD((Xjk + 1 * ldX));
      tz = SCALAR_LOAD((Gik + 0 * ldG));
      tw = SCALAR_LOAD((Gjk + 1 * ldG));
      t0 = SCALAR_LOAD((temp + 0 * blockDim.x + threadIdx.x));
      t0 = SCALAR_MUL(t0, const_value_w);
      tz = SCALAR_FMA(ty, t0, tz);
      tw = SCALAR_FMA(tx, t0, tw);
      SCALAR_STORE((Gik + 0 * ldG), tz);
      SCALAR_STORE((Gjk + 1 * ldG), tw);
      tx = SCALAR_LOAD((Xik + 1 * ldX));
      ty = SCALAR_LOAD((Xjk + 1 * ldX));
      tz = SCALAR_LOAD((Gik + 1 * ldG));
      tw = SCALAR_LOAD((Gjk + 1 * ldG));
      t1 = SCALAR_LOAD((temp + 1 * blockDim.x + threadIdx.x));
      t1 = SCALAR_MUL(t1, const_value_w);
      tz = SCALAR_FMA(ty, t1, tz);
      tw = SCALAR_FMA(tx, t1, tw);
      SCALAR_STORE((Gik + 1 * ldG), tz);
      SCALAR_STORE((Gjk + 1 * ldG), tw);
      tx = SCALAR_LOAD((Xik + 2 * ldX));
      ty = SCALAR_LOAD((Xjk + 1 * ldX));
      tz = SCALAR_LOAD((Gik + 2 * ldG));
      tw = SCALAR_LOAD((Gjk + 1 * ldG));
      t2 = SCALAR_LOAD((temp + 2 * blockDim.x + threadIdx.x));
      t2 = SCALAR_MUL(t2, const_value_w);
      tz = SCALAR_FMA(ty, t2, tz);
      tw = SCALAR_FMA(tx, t2, tw);
      SCALAR_STORE((Gik + 2 * ldG), tz);
      SCALAR_STORE((Gjk + 1 * ldG), tw);
      X_ABp = 1.0; comb_m_i = 1.0;
      Y_ABp = 1.0; comb_n_j = 1.0;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
      const_value_w = SCALAR_MUL(const_value_v, const_value);
      tx = SCALAR_LOAD((Xik + 0 * ldX));
      ty = SCALAR_LOAD((Xjk + 2 * ldX));
      tz = SCALAR_LOAD((Gik + 0 * ldG));
      tw = SCALAR_LOAD((Gjk + 2 * ldG));
      t0 = SCALAR_LOAD((temp + 5 * blockDim.x + threadIdx.x));
      t0 = SCALAR_MUL(t0, const_value_w);
      tz = SCALAR_FMA(ty, t0, tz);
      tw = SCALAR_FMA(tx, t0, tw);
      SCALAR_STORE((Gik + 0 * ldG), tz);
      SCALAR_STORE((Gjk + 2 * ldG), tw);
      tx = SCALAR_LOAD((Xik + 1 * ldX));
      ty = SCALAR_LOAD((Xjk + 2 * ldX));
      tz = SCALAR_LOAD((Gik + 1 * ldG));
      tw = SCALAR_LOAD((Gjk + 2 * ldG));
      t1 = SCALAR_LOAD((temp + 7 * blockDim.x + threadIdx.x));
      t1 = SCALAR_MUL(t1, const_value_w);
      tz = SCALAR_FMA(ty, t1, tz);
      tw = SCALAR_FMA(tx, t1, tw);
      SCALAR_STORE((Gik + 1 * ldG), tz);
      SCALAR_STORE((Gjk + 2 * ldG), tw);
      tx = SCALAR_LOAD((Xik + 2 * ldX));
      ty = SCALAR_LOAD((Xjk + 2 * ldX));
      tz = SCALAR_LOAD((Gik + 2 * ldG));
      tw = SCALAR_LOAD((Gjk + 2 * ldG));
      t2 = SCALAR_LOAD((temp + 8 * blockDim.x + threadIdx.x));
      t2 = SCALAR_MUL(t2, const_value_w);
      tz = SCALAR_FMA(ty, t2, tz);
      tw = SCALAR_FMA(tx, t2, tw);
      SCALAR_STORE((Gik + 2 * ldG), tz);
      SCALAR_STORE((Gjk + 2 * ldG), tw);
      Z_ABp = SCALAR_MUL(Z_ABp, Z_AB); comb_p_k = SCALAR_MUL(comb_p_k * 1, SCALAR_RECIPROCAL(1));
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
      const_value_w = SCALAR_MUL(const_value_v, const_value);
      tx = SCALAR_LOAD((Xik + 0 * ldX));
      ty = SCALAR_LOAD((Xjk + 2 * ldX));
      tz = SCALAR_LOAD((Gik + 0 * ldG));
      tw = SCALAR_LOAD((Gjk + 2 * ldG));
      t0 = SCALAR_LOAD((temp + 0 * blockDim.x + threadIdx.x));
      t0 = SCALAR_MUL(t0, const_value_w);
      tz = SCALAR_FMA(ty, t0, tz);
      tw = SCALAR_FMA(tx, t0, tw);
      SCALAR_STORE((Gik + 0 * ldG), tz);
      SCALAR_STORE((Gjk + 2 * ldG), tw);
      tx = SCALAR_LOAD((Xik + 1 * ldX));
      ty = SCALAR_LOAD((Xjk + 2 * ldX));
      tz = SCALAR_LOAD((Gik + 1 * ldG));
      tw = SCALAR_LOAD((Gjk + 2 * ldG));
      t1 = SCALAR_LOAD((temp + 1 * blockDim.x + threadIdx.x));
      t1 = SCALAR_MUL(t1, const_value_w);
      tz = SCALAR_FMA(ty, t1, tz);
      tw = SCALAR_FMA(tx, t1, tw);
      SCALAR_STORE((Gik + 1 * ldG), tz);
      SCALAR_STORE((Gjk + 2 * ldG), tw);
      tx = SCALAR_LOAD((Xik + 2 * ldX));
      ty = SCALAR_LOAD((Xjk + 2 * ldX));
      tz = SCALAR_LOAD((Gik + 2 * ldG));
      tw = SCALAR_LOAD((Gjk + 2 * ldG));
      t2 = SCALAR_LOAD((temp + 2 * blockDim.x + threadIdx.x));
      t2 = SCALAR_MUL(t2, const_value_w);
      tz = SCALAR_FMA(ty, t2, tz);
      tw = SCALAR_FMA(tx, t2, tw);
      SCALAR_STORE((Gik + 2 * ldG), tz);
      SCALAR_STORE((Gjk + 2 * ldG), tw);
   }
}
