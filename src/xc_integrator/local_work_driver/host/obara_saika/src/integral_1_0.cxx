#include <math.h>
#include "boys_computation.h"
#include "integral_data_types.h"

#define PI 3.14159265358979323846

void integral_1_0(int npts,
                  shells shellA,
                  shells shellB,
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
   double temp[3];

   for(int point_idx = 0; point_idx < npts; ++point_idx) {
      point C = *(_points + point_idx);

      double xA = shellA.origin.x;
      double yA = shellA.origin.y;
      double zA = shellA.origin.z;

      double xB = shellB.origin.x;
      double yB = shellB.origin.y;
      double zB = shellB.origin.z;

      double X_AB = (xA - xB);
      double Y_AB = (yA - yB);
      double Z_AB = (zA - zB);

      double beta_in = 0.0;
      for(int i = 0; i < shellA.m; ++i) {
         for(int j = 0; j < shellB.m; ++j) {
            double aA = shellA.coeff[i].alpha;
            double cA = shellA.coeff[i].coeff;

            double aB = shellB.coeff[j].alpha;
            double cB = shellB.coeff[j].coeff;

            double RHO = aA + aB;
            double RHO_INV = 1.0 / RHO;

            double xC = C.x;
            double yC = C.y;
            double zC = C.z;

            double xP = (aA * xA + aB * xB) * RHO_INV;
            double yP = (aA * yA + aB * yB) * RHO_INV;
            double zP = (aA * zA + aB * zB) * RHO_INV;

            double X_PA = (xP - xA);
            double Y_PA = (yP - yA);
            double Z_PA = (zP - zA);

            double X_PC = (xP - xC);
            double Y_PC = (yP - yC);
            double Z_PC = (zP - zC);

            double t00, t01, t10;

            double eval = cA * cB * 2 * PI * RHO_INV * exp(-1.0 * (X_AB * X_AB + Y_AB * Y_AB + Z_AB * Z_AB) * aA * aB * RHO_INV);
            double tval = RHO * (X_PC * X_PC + Y_PC * Y_PC + Z_PC * Z_PC);

#ifdef BOYS_REFERENCE
            t00 = eval * boys_reference(0, tval);
            t01 = eval * boys_reference(1, tval);
#elif BOYS_ASYMP
            t00 = eval * boys_asymp(0, tval);
            t01 = eval * boys_asymp(1, tval);
#else
            #error "TYPE NOT DEFINED!"
#endif

            t10 = X_PA * t00 - X_PC * t01;
            *(temp + 0) = beta_in * (*(temp + 0)) + t10;

            t10 = Y_PA * t00 - Y_PC * t01;
            *(temp + 1) = beta_in * (*(temp + 1)) + t10;

            t10 = Z_PA * t00 - Z_PC * t01;
            *(temp + 2) = beta_in * (*(temp + 2)) + t10;

            beta_in = 1.0;
         }
      }

      double *Xik = (Xi + point_idx * stX);
      double *Xjk = (Xj + point_idx * stX);
      double *Gik = (Gi + point_idx * stG);
      double *Gjk = (Gj + point_idx * stG);

      double const_value, X_ABp, Y_ABp, Z_ABp, comb_m_i, comb_n_j, comb_p_k, rcp_i, rcp_j, rcp_k;
      double t0, t1, t2;

      X_ABp = 1.0; comb_m_i = 1.0;
      Y_ABp = 1.0; comb_n_j = 1.0;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 0) * const_value;
      *(Gjk + 0 * ldG) += *(Xik + 0 * ldX) * t0;
      *(Gik + 0 * ldG) += *(Xjk + 0 * ldX) * t0;
      t1 = *(temp + 1) * const_value;
      *(Gjk + 1 * ldG) += *(Xik + 0 * ldX) * t1;
      *(Gik + 0 * ldG) += *(Xjk + 1 * ldX) * t1;
      t2 = *(temp + 2) * const_value;
      *(Gjk + 2 * ldG) += *(Xik + 0 * ldX) * t2;
      *(Gik + 0 * ldG) += *(Xjk + 2 * ldX) * t2;

      *(Gjk + 0 * ldG) *= *(weights + point_idx);
      *(Gjk + 1 * ldG) *= *(weights + point_idx);
      *(Gjk + 2 * ldG) *= *(weights + point_idx);

      *(Gik + 0 * ldG) *= *(weights + point_idx);
   }
}
