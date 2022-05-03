#include <stdio.h>
#include <stdlib.h>
#include "../include/gpu/integral_data_types.hpp"
#include "../include/gpu/obara_saika_integrals.hpp"
#include "integral_0.hu"
#include "integral_1.hu"
#include "integral_2.hu"
#include "integral_0_0.hu"
#include "integral_1_0.hu"
#include "integral_1_1.hu"
#include "integral_2_0.hu"
#include "integral_2_1.hu"
#include "integral_2_2.hu"
namespace XGPU {

  void generate_shell_pair( const shells& A, const shells& B, prim_pair *prim_pairs) {
    // L Values
    const auto xA = A.origin.x;
    const auto yA = A.origin.y;
    const auto zA = A.origin.z;

    const auto xB = B.origin.x;
    const auto yB = B.origin.y;
    const auto zB = B.origin.z;

    double rABx = xA - xB;
    double rABy = yA - yB;
    double rABz = zA - zB;

    const double dAB = rABx*rABx + rABy*rABy + rABz*rABz;

    const int nprim_A = A.m;
    const int nprim_B = B.m;
    for(int i = 0, ij = 0; i < nprim_A; ++i       )
      for(int j = 0        ; j < nprim_B; ++j, ++ij ) {
	auto& pair = prim_pairs[ij];
	const auto alpha_A = A.coeff[i].alpha;
	const auto alpha_B = B.coeff[j].alpha;

	pair.gamma = alpha_A + alpha_B;
	pair.gamma_inv = 1. / pair.gamma;

	pair.P.x = (alpha_A * xA + alpha_B * xB) * pair.gamma_inv;
	pair.P.y = (alpha_A * yA + alpha_B * yB) * pair.gamma_inv;
	pair.P.z = (alpha_A * zA + alpha_B * zB) * pair.gamma_inv;

	pair.PA.x = pair.P.x - xA;
	pair.PA.y = pair.P.y - yA;
	pair.PA.z = pair.P.z - zA;

	pair.PB.x = pair.P.x - xB;
	pair.PB.y = pair.P.y - yB;
	pair.PB.z = pair.P.z - zB;

	pair.K_coeff_prod = 2 * M_PI * pair.gamma_inv * std::exp( - alpha_A * alpha_B * dAB * pair.gamma_inv ) * A.coeff[i].coeff * B.coeff[j].coeff;
      }
  }

  void compute_integral_shell_pair(int is_diag,
				   size_t npts,
				   double *points_x,
				   double *points_y,
				   double *points_z,
				   int lA,
				   int lB,
				   point rA,
				   point rB,
           shell_pair* sp,
				   double *Xi,
				   double *Xj,
				   int ldX,
				   double *Gi,
				   double *Gj,
				   int ldG, 
				   double *weights,
		  double *boys_table,
      cudaStream_t stream) {
    if (is_diag) {
      if(lA == 0) {
	integral_0(npts,
		   points_x,
		   points_y,
		   points_z,
       sp,
		   Xi,
		   ldX,
		   Gi,
		   ldG, 
		   weights, 
		   boys_table,
       stream);
      } else if(lA == 1) {
        integral_1(npts,
		   points_x,
		   points_y,
		   points_z,
       sp,
		   Xi,
		   ldX,
		   Gi,
		   ldG, 
		   weights, 
		   boys_table,
       stream);
      } else if(lA == 2) {
        integral_2(npts,
		   points_x,
		   points_y,
		   points_z,
       sp,
		   Xi,
		   ldX,
		   Gi,
		   ldG, 
		   weights, 
		   boys_table,
       stream);
      } else {
	printf("Type not defined!\n");
      }
    } else {
      if((lA == 0) && (lB == 0)) {
	integral_0_0(npts,
		   points_x,
		   points_y,
		   points_z,
         sp,
		     Xi,
		     Xj,
		     ldX,
		     Gi,
		     Gj,
		     ldG, 
		     weights,
		   boys_table,
       stream);
      } else if((lA == 1) && (lB == 0)) {
	integral_1_0(npts,
		   points_x,
		   points_y,
		   points_z,
         sp,
		     Xi,
		     Xj,
		     ldX,
		     Gi,
		     Gj,
		     ldG, 
		     weights,
		   boys_table,
       stream);
      } else if((lA == 0) && (lB == 1)) {
	integral_1_0(npts,
		   points_x,
		   points_y,
		   points_z,
         sp,
		     Xj,
		     Xi,
		     ldX,
		     Gj,
		     Gi,
		     ldG, 
		     weights, 
		   boys_table,
       stream);
      } else if((lA == 1) && (lB == 1)) {
        integral_1_1(rA.x - rB.x,
		     rA.y - rB.y,
		     rA.z - rB.z,
		     npts,
		   points_x,
		   points_y,
		   points_z,
         sp,
		     Xi,
		     Xj,
		     ldX,
		     Gi,
		     Gj,
		     ldG, 
		     weights,
		   boys_table,
       stream);
      } else if((lA == 2) && (lB == 0)) {
	integral_2_0(npts,
		   points_x,
		   points_y,
		   points_z,
         sp,
		     Xi,
		     Xj,
		     ldX,
		     Gi,
		     Gj,
		     ldG, 
		     weights,
		   boys_table,
       stream);
      } else if((lA == 0) && (lB == 2)) {
	integral_2_0(npts,
		   points_x,
		   points_y,
		   points_z,
         sp,
		     Xj,
		     Xi,
		     ldX,
		     Gj,
		     Gi,
		     ldG, 
		     weights, 
		   boys_table,
       stream);
      } else if((lA == 2) && (lB == 1)) {
	integral_2_1(rA.x - rB.x,
		     rA.y - rB.y,
		     rA.z - rB.z,
		     npts,
		   points_x,
		   points_y,
		   points_z,
         sp,
		     Xi,
		     Xj,
		     ldX,
		     Gi,
		     Gj,
		     ldG, 
		     weights,
		   boys_table,
       stream);
      } else if((lA == 1) && (lB == 2)) {
	integral_2_1(rB.x - rA.x,
		     rB.y - rA.y,
		     rB.z - rA.z,
		     npts,
		   points_x,
		   points_y,
		   points_z,
         sp,
		     Xj,
		     Xi,
		     ldX,
		     Gj,
		     Gi,
		     ldG, 
		     weights, 
		   boys_table,
       stream);
      } else if((lA == 2) && (lB == 2)) {
        integral_2_2(rA.x - rB.x,
		     rA.y - rB.y,
		     rA.z - rB.z,
		     npts,
		   points_x,
		   points_y,
		   points_z,
         sp,
		     Xi,
		     Xj,
		     ldX,
		     Gi,
		     Gj,
		     ldG, 
		     weights,
		   boys_table,
       stream);
      } else {
	printf("Type not defined!\n");
      }
    }
  }
}
