#include <math.h>
#include "boys_computation.h"
#include "integral_data_types.h"

#define PI 3.14159265358979323846

#define MIN(a,b)			\
  ({ __typeof__ (a) _a = (a);	        \
  __typeof__ (b) _b = (b);		\
  _a < _b ? _a : _b; })

void integral_0_0(size_t npts,
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
   double temp[1 * NPTS_LOCAL];

   for(int i = 0; i < 1 * NPTS_LOCAL; ++i) {
      temp[i] = 0.0;
   }

   double X_AB = shpair.rAB.x;
   double Y_AB = shpair.rAB.y;
   double Z_AB = shpair.rAB.z;

   for(size_t p_outer = 0; p_outer < npts; p_outer += NPTS_LOCAL) {
      size_t npts_inner = MIN(NPTS_LOCAL, npts - p_outer);
      point *_point_outer = (_points + p_outer);

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

         for(int p_inner = 0; p_inner < npts_inner; ++p_inner) {
            point C = *(_point_outer + p_inner);

            double xC = C.x;
            double yC = C.y;
            double zC = C.z;

            double X_PC = (xP - xC);
            double Y_PC = (yP - yC);
            double Z_PC = (zP - zC);

            double t00;

            double tval = RHO * (X_PC * X_PC + Y_PC * Y_PC + Z_PC * Z_PC);

            t00 = eval * boys_function(0, tval);
            *(temp + 0 * NPTS_LOCAL + p_inner) += t00;
         }
      }

      for(int p_inner = 0; p_inner < npts_inner; ++p_inner) {
         double *Xik = (Xi + (p_outer + p_inner) * stX);
         double *Xjk = (Xj + (p_outer + p_inner) * stX);
         double *Gik = (Gi + (p_outer + p_inner) * stG);
         double *Gjk = (Gj + (p_outer + p_inner) * stG);

         double const_value, X_ABp, Y_ABp, Z_ABp, comb_m_i, comb_n_j, comb_p_k, rcp_i, rcp_j, rcp_k;
         double t0;

         X_ABp = 1.0; comb_m_i = 1.0;
         Y_ABp = 1.0; comb_n_j = 1.0;
         Z_ABp = 1.0; comb_p_k = 1.0;
         const_value = *(weights + p_outer + p_inner) * comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
         t0 = *(temp + 0 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 0 * ldG) += *(Xjk + 0 * ldX) * t0;
         *(Gjk + 0 * ldG) += *(Xik + 0 * ldX) * t0;
      }
   }
}
