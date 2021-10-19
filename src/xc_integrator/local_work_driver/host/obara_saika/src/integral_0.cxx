#include <math.h>
#include "boys_computation.h"
#include "integral_data_types.h"

#define PI 3.14159265358979323846

void integral_0(size_t npts,
               shell_pair shpair,
               point *_points,
               double *Xi,
               int stX,
               int ldX,
               double *Gi,
               int stG, 
               int ldG, 
               double *weights) {
   double temp[1];

   for(int i = 0; i < 1; ++i) {
      temp[i] = 0.0;
   }

   for(size_t point_idx = 0; point_idx < npts; ++point_idx) {
      point C = *(_points + point_idx);

      double xA = shpair.rA.x;
      double yA = shpair.rA.y;
      double zA = shpair.rA.z;

      double beta_in = 0.0;
      for( int ij = 0; ij < shpair.nprim_pair; ++ij ) {
            double RHO = shpair.prim_pairs[ij].gamma;
            double RHO_INV = 1.0 / RHO;

            double xC = C.x;
            double yC = C.y;
            double zC = C.z;

            constexpr double X_PA = 0.0;
            constexpr double Y_PA = 0.0;
            constexpr double Z_PA = 0.0;

            double X_PC = (xA - xC);
            double Y_PC = (yA - yC);
            double Z_PC = (zA - zC);

            double t00;

            double eval = shpair.prim_pairs[ij].coeff_prod * 2 * PI * RHO_INV;
            double tval = RHO * (X_PC * X_PC + Y_PC * Y_PC + Z_PC * Z_PC);

            t00 = eval * boys_function(0, tval);
            *(temp + 0) = beta_in * (*(temp + 0)) + t00;

            beta_in = 1.0;
      }

      double *Xik = (Xi + point_idx * stX);
      double *Gik = (Gi + point_idx * stG);

      double t0;

      t0 = *(temp + 0) * (*(weights + point_idx));
      *(Gik + 0 * ldG) += *(Xik + 0 * ldX) * t0;
   }
}
