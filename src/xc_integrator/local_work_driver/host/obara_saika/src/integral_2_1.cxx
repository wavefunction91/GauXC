#include <math.h>
#include "boys_computation.h"
#include "integral_data_types.h"

#define PI 3.14159265358979323846

#define MIN(a,b)			\
  ({ __typeof__ (a) _a = (a);	        \
  __typeof__ (b) _b = (b);		\
  _a < _b ? _a : _b; })

void integral_2_1(size_t npts,
                  shell_pair shpair,
                  point *_points,
                  double *Xi,
                  double *Xj,
                  int stX,
                  int ldX,
                  double *Gi,
                  double *Gj,
                  int stG, 
                  int ldG, 
                  double *weights) {
   double temp[16 * NPTS_LOCAL];
   double FmT [4 * NPTS_LOCAL];
   double Tval[NPTS_LOCAL];

   double X_AB = shpair.rAB.x;
   double Y_AB = shpair.rAB.y;
   double Z_AB = shpair.rAB.z;

   for(size_t p_outer = 0; p_outer < npts; p_outer += NPTS_LOCAL) {
      size_t npts_inner = MIN(NPTS_LOCAL, npts - p_outer);
      point *_point_outer = (_points + p_outer);

      for(int i = 0; i < 16 * NPTS_LOCAL; ++i) temp[i] = 0.0;

      for(int ij = 0; ij < shpair.nprim_pair; ++ij ) {
         double RHO = shpair.prim_pairs[ij].gamma;
         double RHO_INV = 1.0 / RHO;

         double xP = shpair.prim_pairs[ij].P.x;
         double yP = shpair.prim_pairs[ij].P.y;
         double zP = shpair.prim_pairs[ij].P.z;

         double X_PA = shpair.prim_pairs[ij].PA.x;
         double Y_PA = shpair.prim_pairs[ij].PA.y;
         double Z_PA = shpair.prim_pairs[ij].PA.z;

         double eval = shpair.prim_pairs[ij].coeff_prod * shpair.prim_pairs[ij].K;

         #if 1
         // Evaluate T Values
         for(int p_inner = 0; p_inner < npts_inner; ++p_inner) {
            point C = *(_point_outer + p_inner);

            double xC = C.x;
            double yC = C.y;
            double zC = C.z;

            double X_PC = (xP - xC);
            double Y_PC = (yP - yC);
            double Z_PC = (zP - zC);

            Tval[p_inner] = RHO * (X_PC * X_PC + Y_PC * Y_PC + Z_PC * Z_PC);
         }

         // Evaluate Boys function
         GauXC::boys_function<0>(npts_inner, Tval, FmT + 0*NPTS_LOCAL);
         GauXC::boys_function<1>(npts_inner, Tval, FmT + 1*NPTS_LOCAL);
         GauXC::boys_function<2>(npts_inner, Tval, FmT + 2*NPTS_LOCAL);
         GauXC::boys_function<3>(npts_inner, Tval, FmT + 3*NPTS_LOCAL);

         // Evaluate VRR Buffer
         for(int p_inner = 0; p_inner < npts_inner; ++p_inner) {
            point C = *(_point_outer + p_inner);

            double xC = C.x;
            double yC = C.y;
            double zC = C.z;

            double X_PC = (xP - xC);
            double Y_PC = (yP - yC);
            double Z_PC = (zP - zC);

            double t00, t01, t02, t03, t10, t11, t12, t20, t21, t30;

            t00 = eval * FmT[p_inner + 0*NPTS_LOCAL];
            t01 = eval * FmT[p_inner + 1*NPTS_LOCAL];
            t02 = eval * FmT[p_inner + 2*NPTS_LOCAL];
            t03 = eval * FmT[p_inner + 3*NPTS_LOCAL];
            t10 = X_PA * t00 - X_PC * t01;
            t11 = X_PA * t01 - X_PC * t02;
            t12 = X_PA * t02 - X_PC * t03;
            t20 = X_PA * t10 - X_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = X_PA * t11 - X_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            *(temp + 0 * NPTS_LOCAL + p_inner) += t20;
            t30 = X_PA * t20 - X_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            *(temp + 6 * NPTS_LOCAL + p_inner) += t30;
            t30 = Y_PA * t20 - Y_PC * t21;
            *(temp + 7 * NPTS_LOCAL + p_inner) += t30;
            t30 = Z_PA * t20 - Z_PC * t21;
            *(temp + 8 * NPTS_LOCAL + p_inner) += t30;
            t20 = Y_PA * t10 - Y_PC * t11;
            t21 = Y_PA * t11 - Y_PC * t12;
            *(temp + 1 * NPTS_LOCAL + p_inner) += t20;
            t30 = Y_PA * t20 - Y_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            *(temp + 9 * NPTS_LOCAL + p_inner) += t30;
            t30 = Z_PA * t20 - Z_PC * t21;
            *(temp + 10 * NPTS_LOCAL + p_inner) += t30;
            t20 = Z_PA * t10 - Z_PC * t11;
            t21 = Z_PA * t11 - Z_PC * t12;
            *(temp + 2 * NPTS_LOCAL + p_inner) += t20;
            t30 = Z_PA * t20 - Z_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            *(temp + 11 * NPTS_LOCAL + p_inner) += t30;
            t10 = Y_PA * t00 - Y_PC * t01;
            t11 = Y_PA * t01 - Y_PC * t02;
            t12 = Y_PA * t02 - Y_PC * t03;
            t20 = Y_PA * t10 - Y_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = Y_PA * t11 - Y_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            *(temp + 3 * NPTS_LOCAL + p_inner) += t20;
            t30 = Y_PA * t20 - Y_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            *(temp + 12 * NPTS_LOCAL + p_inner) += t30;
            t30 = Z_PA * t20 - Z_PC * t21;
            *(temp + 13 * NPTS_LOCAL + p_inner) += t30;
            t20 = Z_PA * t10 - Z_PC * t11;
            t21 = Z_PA * t11 - Z_PC * t12;
            *(temp + 4 * NPTS_LOCAL + p_inner) += t20;
            t30 = Z_PA * t20 - Z_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            *(temp + 14 * NPTS_LOCAL + p_inner) += t30;
            t10 = Z_PA * t00 - Z_PC * t01;
            t11 = Z_PA * t01 - Z_PC * t02;
            t12 = Z_PA * t02 - Z_PC * t03;
            t20 = Z_PA * t10 - Z_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = Z_PA * t11 - Z_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            *(temp + 5 * NPTS_LOCAL + p_inner) += t20;
            t30 = Z_PA * t20 - Z_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            *(temp + 15 * NPTS_LOCAL + p_inner) += t30;
         }

         #else

         for(int p_inner = 0; p_inner < npts_inner; ++p_inner) {
            point C = *(_point_outer + p_inner);

            double xC = C.x;
            double yC = C.y;
            double zC = C.z;

            double X_PC = (xP - xC);
            double Y_PC = (yP - yC);
            double Z_PC = (zP - zC);

            double t00, t01, t02, t03, t10, t11, t12, t20, t21, t30;

            double tval = RHO * (X_PC * X_PC + Y_PC * Y_PC + Z_PC * Z_PC);

            t00 = eval * FmT[p_inner + 0*NPTS_LOCAL];
            t01 = eval * FmT[p_inner + 1*NPTS_LOCAL];
            t02 = eval * FmT[p_inner + 2*NPTS_LOCAL];
            t03 = eval * FmT[p_inner + 3*NPTS_LOCAL];
            t10 = X_PA * t00 - X_PC * t01;
            t11 = X_PA * t01 - X_PC * t02;
            t12 = X_PA * t02 - X_PC * t03;
            t20 = X_PA * t10 - X_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = X_PA * t11 - X_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            *(temp + 0 * NPTS_LOCAL + p_inner) += t20;
            t30 = X_PA * t20 - X_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            *(temp + 6 * NPTS_LOCAL + p_inner) += t30;
            t30 = Y_PA * t20 - Y_PC * t21;
            *(temp + 7 * NPTS_LOCAL + p_inner) += t30;
            t30 = Z_PA * t20 - Z_PC * t21;
            *(temp + 8 * NPTS_LOCAL + p_inner) += t30;
            t20 = Y_PA * t10 - Y_PC * t11;
            t21 = Y_PA * t11 - Y_PC * t12;
            *(temp + 1 * NPTS_LOCAL + p_inner) += t20;
            t30 = Y_PA * t20 - Y_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            *(temp + 9 * NPTS_LOCAL + p_inner) += t30;
            t30 = Z_PA * t20 - Z_PC * t21;
            *(temp + 10 * NPTS_LOCAL + p_inner) += t30;
            t20 = Z_PA * t10 - Z_PC * t11;
            t21 = Z_PA * t11 - Z_PC * t12;
            *(temp + 2 * NPTS_LOCAL + p_inner) += t20;
            t30 = Z_PA * t20 - Z_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            *(temp + 11 * NPTS_LOCAL + p_inner) += t30;
            t10 = Y_PA * t00 - Y_PC * t01;
            t11 = Y_PA * t01 - Y_PC * t02;
            t12 = Y_PA * t02 - Y_PC * t03;
            t20 = Y_PA * t10 - Y_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = Y_PA * t11 - Y_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            *(temp + 3 * NPTS_LOCAL + p_inner) += t20;
            t30 = Y_PA * t20 - Y_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            *(temp + 12 * NPTS_LOCAL + p_inner) += t30;
            t30 = Z_PA * t20 - Z_PC * t21;
            *(temp + 13 * NPTS_LOCAL + p_inner) += t30;
            t20 = Z_PA * t10 - Z_PC * t11;
            t21 = Z_PA * t11 - Z_PC * t12;
            *(temp + 4 * NPTS_LOCAL + p_inner) += t20;
            t30 = Z_PA * t20 - Z_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            *(temp + 14 * NPTS_LOCAL + p_inner) += t30;
            t10 = Z_PA * t00 - Z_PC * t01;
            t11 = Z_PA * t01 - Z_PC * t02;
            t12 = Z_PA * t02 - Z_PC * t03;
            t20 = Z_PA * t10 - Z_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = Z_PA * t11 - Z_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            *(temp + 5 * NPTS_LOCAL + p_inner) += t20;
            t30 = Z_PA * t20 - Z_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            *(temp + 15 * NPTS_LOCAL + p_inner) += t30;
         }

         #endif
      }

      for(int p_inner = 0; p_inner < npts_inner; ++p_inner) {
         double *Xik = (Xi + (p_outer + p_inner) * stX);
         double *Xjk = (Xj + (p_outer + p_inner) * stX);
         double *Gik = (Gi + (p_outer + p_inner) * stG);
         double *Gjk = (Gj + (p_outer + p_inner) * stG);

         double const_value, X_ABp, Y_ABp, Z_ABp, comb_m_i, comb_n_j, comb_p_k, rcp_i, rcp_j, rcp_k;
         double t0, t1, t2, t3, t4, t5;

         X_ABp = 1.0; comb_m_i = 1.0;
         Y_ABp = 1.0; comb_n_j = 1.0;
         Z_ABp = 1.0; comb_p_k = 1.0;
         const_value = *(weights + p_outer + p_inner) * comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
         t0 = *(temp + 6 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 0 * ldG) += *(Xjk + 0 * ldX) * t0;
         *(Gjk + 0 * ldG) += *(Xik + 0 * ldX) * t0;
         t1 = *(temp + 7 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 1 * ldG) += *(Xjk + 0 * ldX) * t1;
         *(Gjk + 0 * ldG) += *(Xik + 1 * ldX) * t1;
         t2 = *(temp + 8 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 2 * ldG) += *(Xjk + 0 * ldX) * t2;
         *(Gjk + 0 * ldG) += *(Xik + 2 * ldX) * t2;
         t3 = *(temp + 9 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 3 * ldG) += *(Xjk + 0 * ldX) * t3;
         *(Gjk + 0 * ldG) += *(Xik + 3 * ldX) * t3;
         t4 = *(temp + 10 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 4 * ldG) += *(Xjk + 0 * ldX) * t4;
         *(Gjk + 0 * ldG) += *(Xik + 4 * ldX) * t4;
         t5 = *(temp + 11 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 5 * ldG) += *(Xjk + 0 * ldX) * t5;
         *(Gjk + 0 * ldG) += *(Xik + 5 * ldX) * t5;
         X_ABp *= X_AB; rcp_i = 1.0 / (1.0 * 1); comb_m_i = (comb_m_i * 1) * rcp_i;
         Y_ABp = 1.0; comb_n_j = 1.0;
         Z_ABp = 1.0; comb_p_k = 1.0;
         const_value = *(weights + p_outer + p_inner) * comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
         t0 = *(temp + 0 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 0 * ldG) += *(Xjk + 0 * ldX) * t0;
         *(Gjk + 0 * ldG) += *(Xik + 0 * ldX) * t0;
         t1 = *(temp + 1 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 1 * ldG) += *(Xjk + 0 * ldX) * t1;
         *(Gjk + 0 * ldG) += *(Xik + 1 * ldX) * t1;
         t2 = *(temp + 2 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 2 * ldG) += *(Xjk + 0 * ldX) * t2;
         *(Gjk + 0 * ldG) += *(Xik + 2 * ldX) * t2;
         t3 = *(temp + 3 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 3 * ldG) += *(Xjk + 0 * ldX) * t3;
         *(Gjk + 0 * ldG) += *(Xik + 3 * ldX) * t3;
         t4 = *(temp + 4 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 4 * ldG) += *(Xjk + 0 * ldX) * t4;
         *(Gjk + 0 * ldG) += *(Xik + 4 * ldX) * t4;
         t5 = *(temp + 5 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 5 * ldG) += *(Xjk + 0 * ldX) * t5;
         *(Gjk + 0 * ldG) += *(Xik + 5 * ldX) * t5;
         X_ABp = 1.0; comb_m_i = 1.0;
         Y_ABp = 1.0; comb_n_j = 1.0;
         Z_ABp = 1.0; comb_p_k = 1.0;
         const_value = *(weights + p_outer + p_inner) * comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
         t0 = *(temp + 7 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 0 * ldG) += *(Xjk + 1 * ldX) * t0;
         *(Gjk + 1 * ldG) += *(Xik + 0 * ldX) * t0;
         t1 = *(temp + 9 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 1 * ldG) += *(Xjk + 1 * ldX) * t1;
         *(Gjk + 1 * ldG) += *(Xik + 1 * ldX) * t1;
         t2 = *(temp + 10 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 2 * ldG) += *(Xjk + 1 * ldX) * t2;
         *(Gjk + 1 * ldG) += *(Xik + 2 * ldX) * t2;
         t3 = *(temp + 12 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 3 * ldG) += *(Xjk + 1 * ldX) * t3;
         *(Gjk + 1 * ldG) += *(Xik + 3 * ldX) * t3;
         t4 = *(temp + 13 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 4 * ldG) += *(Xjk + 1 * ldX) * t4;
         *(Gjk + 1 * ldG) += *(Xik + 4 * ldX) * t4;
         t5 = *(temp + 14 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 5 * ldG) += *(Xjk + 1 * ldX) * t5;
         *(Gjk + 1 * ldG) += *(Xik + 5 * ldX) * t5;
         Y_ABp *= Y_AB; rcp_j = 1.0 / (1.0 * 1); comb_n_j = (comb_n_j * 1) * rcp_j;
         Z_ABp = 1.0; comb_p_k = 1.0;
         const_value = *(weights + p_outer + p_inner) * comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
         t0 = *(temp + 0 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 0 * ldG) += *(Xjk + 1 * ldX) * t0;
         *(Gjk + 1 * ldG) += *(Xik + 0 * ldX) * t0;
         t1 = *(temp + 1 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 1 * ldG) += *(Xjk + 1 * ldX) * t1;
         *(Gjk + 1 * ldG) += *(Xik + 1 * ldX) * t1;
         t2 = *(temp + 2 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 2 * ldG) += *(Xjk + 1 * ldX) * t2;
         *(Gjk + 1 * ldG) += *(Xik + 2 * ldX) * t2;
         t3 = *(temp + 3 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 3 * ldG) += *(Xjk + 1 * ldX) * t3;
         *(Gjk + 1 * ldG) += *(Xik + 3 * ldX) * t3;
         t4 = *(temp + 4 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 4 * ldG) += *(Xjk + 1 * ldX) * t4;
         *(Gjk + 1 * ldG) += *(Xik + 4 * ldX) * t4;
         t5 = *(temp + 5 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 5 * ldG) += *(Xjk + 1 * ldX) * t5;
         *(Gjk + 1 * ldG) += *(Xik + 5 * ldX) * t5;
         X_ABp = 1.0; comb_m_i = 1.0;
         Y_ABp = 1.0; comb_n_j = 1.0;
         Z_ABp = 1.0; comb_p_k = 1.0;
         const_value = *(weights + p_outer + p_inner) * comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
         t0 = *(temp + 8 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 0 * ldG) += *(Xjk + 2 * ldX) * t0;
         *(Gjk + 2 * ldG) += *(Xik + 0 * ldX) * t0;
         t1 = *(temp + 10 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 1 * ldG) += *(Xjk + 2 * ldX) * t1;
         *(Gjk + 2 * ldG) += *(Xik + 1 * ldX) * t1;
         t2 = *(temp + 11 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 2 * ldG) += *(Xjk + 2 * ldX) * t2;
         *(Gjk + 2 * ldG) += *(Xik + 2 * ldX) * t2;
         t3 = *(temp + 13 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 3 * ldG) += *(Xjk + 2 * ldX) * t3;
         *(Gjk + 2 * ldG) += *(Xik + 3 * ldX) * t3;
         t4 = *(temp + 14 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 4 * ldG) += *(Xjk + 2 * ldX) * t4;
         *(Gjk + 2 * ldG) += *(Xik + 4 * ldX) * t4;
         t5 = *(temp + 15 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 5 * ldG) += *(Xjk + 2 * ldX) * t5;
         *(Gjk + 2 * ldG) += *(Xik + 5 * ldX) * t5;
         Z_ABp *= Z_AB; rcp_k = 1.0 / (1.0 * 1); comb_p_k = (comb_p_k * 1) * rcp_k;
         const_value = *(weights + p_outer + p_inner) * comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
         t0 = *(temp + 0 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 0 * ldG) += *(Xjk + 2 * ldX) * t0;
         *(Gjk + 2 * ldG) += *(Xik + 0 * ldX) * t0;
         t1 = *(temp + 1 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 1 * ldG) += *(Xjk + 2 * ldX) * t1;
         *(Gjk + 2 * ldG) += *(Xik + 1 * ldX) * t1;
         t2 = *(temp + 2 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 2 * ldG) += *(Xjk + 2 * ldX) * t2;
         *(Gjk + 2 * ldG) += *(Xik + 2 * ldX) * t2;
         t3 = *(temp + 3 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 3 * ldG) += *(Xjk + 2 * ldX) * t3;
         *(Gjk + 2 * ldG) += *(Xik + 3 * ldX) * t3;
         t4 = *(temp + 4 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 4 * ldG) += *(Xjk + 2 * ldX) * t4;
         *(Gjk + 2 * ldG) += *(Xik + 4 * ldX) * t4;
         t5 = *(temp + 5 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 5 * ldG) += *(Xjk + 2 * ldX) * t5;
         *(Gjk + 2 * ldG) += *(Xik + 5 * ldX) * t5;
      }
   }
}
