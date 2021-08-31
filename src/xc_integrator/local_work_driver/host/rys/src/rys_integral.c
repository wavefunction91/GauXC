#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>

#include "rys_rw.h"
#include "rys_integral.h"

#define Lx 8
#define Ly 8

#define Vx (Lx + 1) * (Lx + 2) / 2
#define Vy (Ly + 1) * (Ly + 2) / 2

#define R_MAX (Lx + 1)

#define PB 16

#define PI 3.14159265358979323846

#define MIN(a,b)                                \
  ({ __typeof__ (a) _a = (a);                   \
    __typeof__ (b) _b = (b);                    \
    _a < _b ? _a : _b; })

// codelets
inline void __attribute__((always_inline)) compute_00(double beta, double *int_array, double *wgh) {
  *(int_array + 0) = (*(int_array + 0)) * beta + *(wgh + 0);
}

inline void __attribute__((always_inline)) compute_10_01(double xPX, double yPX, double zPX, double xPC, double yPC, double zPC, double beta, double *int_array, double *rts, double *wgh) {
  double rt, Cx0, Cy0, Cz0, Cx1, Cy1, Cz1;
  
  rt = *(rts + 0);
  Cx0 = (xPX - xPC * rt);
  Cy0 = (yPX - yPC * rt);
  Cz0 = (zPX - zPC * rt);

  rt = *(rts + 1);
  Cx1 = (xPX - xPC * rt);
  Cy1 = (yPX - yPC * rt);
  Cz1 = (zPX - zPC * rt);

  *(int_array + 0) = (*(int_array + 0)) * beta + (*(wgh + 0)) * Cx0 + (*(wgh + 1)) * Cx1;
  *(int_array + 1) = (*(int_array + 1)) * beta + (*(wgh + 0)) * Cy0 + (*(wgh + 1)) * Cy1;
  *(int_array + 2) = (*(int_array + 2)) * beta + (*(wgh + 0)) * Cz0 + (*(wgh + 1)) * Cz1;
}

inline void __attribute__((always_inline)) compute_20_02(double xPX, double yPX, double zPX, double xPC, double yPC, double zPC, double aP_inv, double beta, double *int_array, double *rts, double *wgh) {
  double B0, B1, rt0, rt1, Cx0, Cy0, Cz0, Cx1, Cy1, Cz1, Cx2, Cy2, Cz2, Cx3, Cy3, Cz3;
  
  rt0 = *(rts + 0);
  Cx0 = (xPX - xPC * rt0);
  Cy0 = (yPX - yPC * rt0);
  Cz0 = (zPX - zPC * rt0);

  rt1 = *(rts + 1);
  Cx1 = (xPX - xPC * rt1);
  Cy1 = (yPX - yPC * rt1);
  Cz1 = (zPX - zPC * rt1);

  B0 = (1.0 - rt0) * aP_inv * 0.5;
  B1 = (1.0 - rt1) * aP_inv * 0.5;
	      
  Cx2 = Cx0 * Cx0 + B0;
  Cy2 = Cy0 * Cy0 + B0;
  Cz2 = Cz0 * Cz0 + B0;

  Cx3 = Cx1 * Cx1 + B1;
  Cy3 = Cy1 * Cy1 + B1;
  Cz3 = Cz1 * Cz1 + B1;
	      
  *(int_array + 0) = (*(int_array + 0)) * beta + Cx2 * (*(wgh + 0)) + Cx3 * (*(wgh + 1));
  *(int_array + 1) = (*(int_array + 1)) * beta + Cx0 * Cy0 * (*(wgh + 0)) + Cx1 * Cy1 * (*(wgh + 1));
  *(int_array + 2) = (*(int_array + 2)) * beta + Cx0 * Cz0 * (*(wgh + 0)) + Cx1 * Cz1 * (*(wgh + 1));

  *(int_array + 3) = (*(int_array + 3)) * beta + Cy2 * (*(wgh + 0)) + Cy3 * (*(wgh + 1));
  *(int_array + 4) = (*(int_array + 4)) * beta + Cy0 * Cz0 * (*(wgh + 0)) + Cy1 * Cz1 * (*(wgh + 1));
  *(int_array + 5) = (*(int_array + 5)) * beta + Cz2 * (*(wgh + 0)) + Cz3 * (*(wgh + 1));
}

inline void __attribute__((always_inline)) compute_11(double xAB, double yAB, double zAB, double xPX, double yPX, double zPX, double xPC, double yPC, double zPC, double aP_inv, double beta, double *int_array, double *rts, double *wgh) {
  double B0, B1, rt0, rt1, Cx0, Cy0, Cz0, Cx1, Cy1, Cz1, Cx2, Cy2, Cz2, Cx3, Cy3, Cz3;
  
  rt0 = *(rts + 0);
  Cx0 = (xPX - xPC * rt0);
  Cy0 = (yPX - yPC * rt0);
  Cz0 = (zPX - zPC * rt0);

  rt1 = *(rts + 1);
  Cx1 = (xPX - xPC * rt1);
  Cy1 = (yPX - yPC * rt1);
  Cz1 = (zPX - zPC * rt1);

  B0 = (1.0 - rt0) * aP_inv * 0.5;
  B1 = (1.0 - rt1) * aP_inv * 0.5;
	      
  Cx2 = Cx0 * Cx0 + B0;
  Cy2 = Cy0 * Cy0 + B0;
  Cz2 = Cz0 * Cz0 + B0;

  Cx3 = Cx1 * Cx1 + B1;
  Cy3 = Cy1 * Cy1 + B1;
  Cz3 = Cz1 * Cz1 + B1;
	      
  *(int_array + 0) = (*(int_array + 0)) * beta + (Cx2 + xAB * Cx0) * (*(wgh + 0)) + (Cx3 + xAB * Cx1) * (*(wgh + 1));
  *(int_array + 1) = (*(int_array + 1)) * beta + Cy0 * (Cx0 + xAB) * (*(wgh + 0)) + Cy1 * (Cx1 + xAB) * (*(wgh + 1));
  *(int_array + 2) = (*(int_array + 2)) * beta + Cz0 * (Cx0 + xAB) * (*(wgh + 0)) + Cz1 * (Cx1 + xAB) * (*(wgh + 1));

  *(int_array + 3) = (*(int_array + 3)) * beta + Cx0 * (Cy0 + yAB) * (*(wgh + 0)) + Cx1 * (Cy1 + yAB) * (*(wgh + 1));
  *(int_array + 4) = (*(int_array + 4)) * beta + (Cy2 + yAB * Cy0) * (*(wgh + 0)) + (Cy3 + yAB * Cy1) * (*(wgh + 1));
  *(int_array + 5) = (*(int_array + 5)) * beta + Cz0 * (Cy0 + yAB) * (*(wgh + 0)) + Cz1 * (Cy1 + yAB) * (*(wgh + 1));

  *(int_array + 6) = (*(int_array + 6)) * beta + Cx0 * (Cz0 + zAB) * (*(wgh + 0)) + Cx1 * (Cz1 + zAB) * (*(wgh + 1));
  *(int_array + 7) = (*(int_array + 7)) * beta + Cy0 * (Cz0 + zAB) * (*(wgh + 0)) + Cy1 * (Cz1 + zAB) * (*(wgh + 1));
  *(int_array + 8) = (*(int_array + 8)) * beta + (Cz2 + zAB * Cz0) * (*(wgh + 0)) + (Cz3 + zAB * Cz1) * (*(wgh + 1));
}

// nr roots > 2
inline void __attribute__((always_inline)) compute_vrr3(int nr_roots, int l, int lA, int llA, int lB, int llB, double xPX, double yPX, double zPX, double xPC, double yPC, double zPC, double aP_inv, double * rts, double *vrr_array, double *hrr_array) {
  double *roots = (rts + 0);
  double *vrr = (vrr_array + 0);
  for(int r = 0; r < nr_roots; ++r) {
    double *hrr = (hrr_array + l * r);
    
    double rt = *(roots + 0);
    
    double B = (1.0 - rt) * aP_inv * 0.5;

    double Cx = (xPX - xPC * rt);
    double Cy = (yPX - yPC * rt);
    double Cz = (zPX - zPC * rt);

    double v0x = B;
    double v0y = B;
    double v0z = B;

    double v1x = Cx;
    double v1y = Cy;
    double v1z = Cz;
    
    *(vrr + 0) = 1.0;
    *(vrr + 1) = 1.0;
    *(vrr + 2) = 1.0;
    
    *(vrr + 3) = Cx;
    *(vrr + 4) = Cy;
    *(vrr + 5) = Cz;

    *(hrr + 0 * llB + 0) = 1.0;
    *(hrr + 0 * llB + 1) = 1.0;
    *(hrr + 0 * llB + 2) = 1.0;

    *(hrr + 1 * llB + 0) = Cx;
    *(hrr + 1 * llB + 1) = Cy;
    *(hrr + 1 * llB + 2) = Cz;
    
    double *vrri = (vrr + 6);
    double *hrri = (hrr + 2 * llB);
    for(int i = 1; i <= lB; ++i) {
      double v2x = Cx * v1x + i * v0x;
      double v2y = Cy * v1y + i * v0y;
      double v2z = Cz * v1z + i * v0z;
	
      *(vrri + 0) = v2x;
      *(vrri + 1) = v2y;
      *(vrri + 2) = v2z;

      *(hrri + 0) = v2x;
      *(hrri + 1) = v2y;
      *(hrri + 2) = v2z;

      v0x = v1x * B;
      v0y = v1y * B;
      v0z = v1z * B;
	
      v1x = v2x;
      v1y = v2y;
      v1z = v2z;

      vrri += 3;
      hrri += llB;
    }

    for(int i = lB + 1; i <= lA + lB; ++i) {
      double v2x = Cx * v1x + i * v0x;
      double v2y = Cy * v1y + i * v0y;
      double v2z = Cz * v1z + i * v0z;
	
      *(vrri + 0) = v2x;
      *(vrri + 1) = v2y;
      *(vrri + 2) = v2z;

      v0x = v1x * B;
      v0y = v1y * B;
      v0z = v1z * B;
	
      v1x = v2x;
      v1y = v2y;
      v1z = v2z;

      vrri += 3;
    }
    
    vrr += 3 * (lA + lB + 1);
    roots ++;
  }
}

inline void __attribute__((always_inline)) compute_hrr3(int nr_roots, int l, int lA, int llA, int lB, int llB, double xAB, double yAB, double zAB, double *vrr_array, double *hrr_array) {
  for(int j = 1; j <= lA; ++j) {
    double *hrrj = (hrr_array + llA * j);
    
    for(int r = 0; r < nr_roots; ++r) {
      double *hrrr = (hrrj + l * r);
      double *vrrr = (vrr_array + 3 * (lA + lB + 1) * r);
      
      double v0x = *(vrrr + 0);
      double v0y = *(vrrr + 1);
      double v0z = *(vrrr + 2);
      
      for(int i = 0; i <= lB; ++i) {
	double v1x = *(vrrr + 3);
	double v1y = *(vrrr + 4);
	double v1z = *(vrrr + 5);
	
	v0x = v1x + xAB * v0x;
	v0y = v1y + yAB * v0y;
	v0z = v1z + zAB * v0z;

	*(hrrr + 0) = v0x;
	*(hrrr + 1) = v0y;
	*(hrrr + 2) = v0z;

	*(vrrr + 0) = v0x;
	*(vrrr + 1) = v0y;
	*(vrrr + 2) = v0z;

	v0x = v1x;
	v0y = v1y;
	v0z = v1z;

	vrrr += 3;
	hrrr += llB;
      }

      for(int i = lB + 1; i <= lA + lB - j; ++i) {
	double v1x = *(vrrr + 3);
	double v1y = *(vrrr + 4);
	double v1z = *(vrrr + 5);
	
	v0x = v1x + xAB * v0x;
	v0y = v1y + yAB * v0y;
	v0z = v1z + zAB * v0z;

	*(vrrr + 0) = v0x;
	*(vrrr + 1) = v0y;
	*(vrrr + 2) = v0z;

	v0x = v1x;
	v0y = v1y;
	v0z = v1z;

	vrrr += 3;
	hrrr += llB;
      }
    }
  }
}

inline int __attribute__((always_inline)) index_calculation(int i, int j, int L) {
  return (L - i) * (L - i + 1) / 2 + j;
}

inline void __attribute__((always_inline)) compute_reduction(int nr_roots, int lA, int lB, double *weights, double *hrr_array, double *result, double beta) {
  int offsetB = (lB + 1) * (lB + 2) / 2;

  for(int ia = 0; ia <= lA; ++ia) {
    for(int ja = 0; ja <= (lA - ia); ++ja) {
      int ka = lA - ia - ja;
      int ija = index_calculation(ia, ka, lA);

      double *hrria = (hrr_array + 3 * (lB + 1) * nr_roots * ia);
      double *hrrja = (hrr_array + 3 * (lB + 1) * nr_roots * ja);
      double *hrrka = (hrr_array + 3 * (lB + 1) * nr_roots * ka);
      
      for(int ib = 0; ib <= lB; ++ib) {
	for(int jb = 0; jb <= (lB - ib); ++jb) {
	  int kb = lB - ib - jb;
	  int ijb = index_calculation(ib, kb, lB);

	  double *hrrib = (hrria + 3 * ib);
	  double *hrrjb = (hrrja + 3 * jb);
	  double *hrrkb = (hrrka + 3 * kb);
	  double *wghs = (weights + 0);
	  
	  double value = 0.0;

	  for(int r = 0; r < nr_roots; ++r) {
	    double ix = *(hrrib + 0);
	    double iy = *(hrrjb + 1);
	    double iz = *(hrrkb + 2);
	    double w  = *(wghs + 0);

	    value += (ix * iy * iz * w);

	    hrrib += 3 * (lB + 1);
	    hrrjb += 3 * (lB + 1);
	    hrrkb += 3 * (lB + 1);
	    wghs ++;
	  }

	  *(result + offsetB * ija + ijb) = (*(result + offsetB * ija + ijb)) * beta + value;
	}
      }
    }
  }
}

void compute_integral(int n, shells *shell_list, int m, point *points, double *matrix) {
  double *rts = (double*) malloc(PB * R_MAX * sizeof(double));
  double *wgh = (double*) malloc(PB * R_MAX * sizeof(double));

  double *int_array = (double*) malloc(PB * Vx * Vy * sizeof(double));
  double *vrr_array = (double*) malloc(3 * (Lx + Ly + 1) * R_MAX * sizeof(double));
  double *hrr_array = (double*) malloc(3 * (Lx + 1) * (Ly + 1) * R_MAX * sizeof(double));

  int nn = 0;
  for(int i = 0; i < n; ++i) {
    int L = shell_list[i].L;    
    nn += ((L + 1) * (L + 2) / 2);
  }
  
  int offset_ii = 0;
  for(int ii = 0; ii < n; ++ii) {
    shells shell0 = shell_list[ii];

    int offset_jj = 0;
    for(int jj = 0; jj < n; ++jj) {
      shells shell1 = shell_list[jj];

      for(int p = 0; p < m; p += PB) {
	int pp = MIN(m - p, PB);
	point *ppoints = (points + p);
      
	// values
	double xA = shell0.origin.x;
	double yA = shell0.origin.y;
	double zA = shell0.origin.z;
	int lA = shell0.L;
	    
	double xB = shell1.origin.x;
	double yB = shell1.origin.y;
	double zB = shell1.origin.z;
	int lB = shell1.L;

	// nr of roots
	int nr_roots = ((int) ceil((lA + lB) / 2.0)) + 1;

	double xAB = (lB < lA) ? (xA - xB) : (xB - xA);
	double yAB = (lB < lA) ? (yA - yB) : (yB - yA);
	double zAB = (lB < lA) ? (zA - zB) : (zB - zA);
	
	double beta = 0.0;
	for(int i = 0; i < shell0.m; ++i) {
	  for(int j = 0; j < shell1.m; ++j) {
	    
	    double aA = shell0.coeff[i].alpha;
	    double cA = shell0.coeff[i].coeff;  

	    double aB = shell1.coeff[j].alpha;
	    double cB = shell1.coeff[j].coeff;

	    double aP = aA + aB;
	    double aP_inv = 1.0 / aP;
  
	    double xP = (aA * xA + aB * xB) * aP_inv;
	    double yP = (aA * yA + aB * yB) * aP_inv;
	    double zP = (aA * zA + aB * zB) * aP_inv;
  
	    double xPX = (lB < lA) ? (xP - xA) : (xP - xB);
	    double yPX = (lB < lA) ? (yP - yA) : (yP - yB);
	    double zPX = (lB < lA) ? (zP - zA) : (zP - zB);

	    double tval[PB];
	    double xPC[PB];
	    double yPC[PB];
	    double zPC[PB];
	    
	    double eval = exp(-1.0 * (xAB * xAB + yAB * yAB + zAB * zAB) * aA * aB * aP_inv);

	    for(int pb = 0; pb < pp; ++pb) {
	      point C = *(ppoints + pb);
	      
	      double xC = (xP - C.x);
	      double yC = (yP - C.y);
	      double zC = (zP - C.z);

	      xPC[pb] = xC;
	      yPC[pb] = yC;
	      zPC[pb] = zC;
	      
	      tval[pb] = aP * (xC * xC + yC * yC + zC * zC);
	    }

	    for(int pb = 0; pb < pp * nr_roots; ++pb) {
	      *(rts + pb) = 0.0;
	      *(wgh + pb) = 2 * PI * aP_inv * eval * cA * cB;
	    }
  
	    rys_rw(pp, nr_roots, tval, rts, wgh);  

	    int lX = (lB < lA) ? lB : lA;
	    int lY = (lB < lA) ? lA : lB;
	    int llX = (lB < lA) ? 3 : 3 * (lB + 1) * nr_roots;
	    int llY = (lB < lA) ? 3 * (lB + 1) * nr_roots : 3;

	    if((lA == 0) && (lB == 0)) {
	      for(int pb = 0; pb < pp; ++pb) {
		compute_00(beta, (int_array + 1 * pb), (wgh + 1 * pb));
	      }
	    } else if(((lA == 1) && (lB == 0)) || ((lA == 0) && (lB == 1))) {
	      for(int pb = 0; pb < pp; ++pb) {
		double xC = xPC[pb];
		double yC = yPC[pb];
		double zC = zPC[pb];
		
		compute_10_01(xPX, yPX, zPX, xC, yC, zC, beta, (int_array + 3 * pb), (rts + 2 * pb), (wgh + 2 * pb));
	      }
	    } else if((lA == 1) && (lB == 1)) {
	      for(int pb = 0; pb < pp; ++pb) {
		double xC = xPC[pb];
		double yC = yPC[pb];
		double zC = zPC[pb];

		compute_11(xAB, yAB, zAB, xPX, yPX, zPX, xC, yC, zC, aP_inv, beta, (int_array + 9 * pb), (rts + 2 * pb), (wgh + 2 * pb));
	      }
	    } else if(((lA == 2) && (lB == 0)) || ((lA == 0) && (lB == 2))) {
	      for(int pb = 0; pb < pp; ++pb) {
		double xC = xPC[pb];
		double yC = yPC[pb];
		double zC = zPC[pb];

		compute_20_02(xPX, yPX, zPX, xC, yC, zC, aP_inv, beta, (int_array + 6 * pb), (rts + 2 * pb), (wgh + 2 * pb));
	      }
	    } else {
	      for(int pb = 0; pb < pp; ++pb) {
		double xC = xPC[pb];
		double yC = yPC[pb];
		double zC = zPC[pb];

		compute_vrr3(nr_roots, 3 * (lB + 1), lX, llX, lY, llY, xPX, yPX, zPX, xC, yC, zC, aP_inv, (rts + nr_roots * pb), (vrr_array + 0), (hrr_array + 0));	    
		compute_hrr3(nr_roots, 3 * (lB + 1), lX, llX, lY, llY, xAB, yAB, zAB, (vrr_array + 0), (hrr_array + 0));
		compute_reduction(nr_roots, lA, lB, (wgh + nr_roots * pb), (hrr_array + 0), (int_array + ((lA + 1) * (lA + 2) / 2) * ((lB + 1) * (lB + 2) / 2) * pb), beta);
	      }
	    }
	    
	    beta = 1.0;
	  }
	}

	for(int pb = 0; pb < pp; ++pb) {
	  for(int i = 0; i < (lA + 1) * (lA + 2) / 2; ++i) {
	    for(int j = 0; j < (lB + 1) * (lB + 2) / 2; ++j) {
	      *(matrix + nn * nn * (p + pb) + nn * (i + offset_ii) + (j + offset_jj)) = *(int_array + ((lA + 1) * (lA + 2) / 2) * ((lB + 1) * (lB + 2) / 2) * pb + ((lB + 1) * (lB + 2) / 2) * i + j);
	    }
	  }
	}
      }
      
      int lB = shell1.L;
      offset_jj += ((lB + 1) * (lB + 2) / 2);
    }

    int lA = shell0.L;
    offset_ii += ((lA + 1) * (lA + 2) / 2);
  }
  
  free(rts);
  free(wgh);
  
  free(int_array);
  free(vrr_array);
  free(hrr_array);
}

void compute_integral_shell_pair( int npts,
				  shells sh0,
				  shells sh1, 
                                  point *points,
				  double *matrix ) {
  double *rts = (double*) malloc(PB * R_MAX * sizeof(double));
  double *wgh = (double*) malloc(PB * R_MAX * sizeof(double));

  double *vrr_array = (double*) malloc(3 * (Lx + Ly + 1) * R_MAX * sizeof(double));
  double *hrr_array = (double*) malloc(3 * (Lx + 1) * (Ly + 1) * R_MAX * sizeof(double));

  // values
  double xA = sh0.origin.x;
  double yA = sh0.origin.y;
  double zA = sh0.origin.z;
  int lA = sh0.L;
	    
  double xB = sh1.origin.x;
  double yB = sh1.origin.y;
  double zB = sh1.origin.z;
  int lB = sh1.L;

  // nr of roots
  int nr_roots = ((int) ceil((lA + lB) / 2.0)) + 1;

  double xAB = (lB < lA) ? (xA - xB) : (xB - xA);
  double yAB = (lB < lA) ? (yA - yB) : (yB - yA);
  double zAB = (lB < lA) ? (zA - zB) : (zB - zA);

  const int shpair_sz =  (lA+1)*(lA+2) * (lB+1)*(lB+2) / 4;
  for(int p = 0; p < npts; p += PB) {
    int pp = MIN(npts - p, PB);
    point *ppoints = (points + p);
    
    double beta = 0.0;
    for(int i = 0; i < sh0.m; ++i) {
      for(int j = 0; j < sh1.m; ++j) {
	    
	double aA = sh0.coeff[i].alpha;
	double cA = sh0.coeff[i].coeff;  

	double aB = sh1.coeff[j].alpha;
	double cB = sh1.coeff[j].coeff;

	double aP = aA + aB;
	double aP_inv = 1.0 / aP;
  
	double xP = (aA * xA + aB * xB) * aP_inv;
	double yP = (aA * yA + aB * yB) * aP_inv;
	double zP = (aA * zA + aB * zB) * aP_inv;
  
	double xPX = (lB < lA) ? (xP - xA) : (xP - xB);
	double yPX = (lB < lA) ? (yP - yA) : (yP - yB);
	double zPX = (lB < lA) ? (zP - zA) : (zP - zB);

	double tval[PB];
	double xPC[PB];
	double yPC[PB];
	double zPC[PB];
	    
	double eval = exp(-1.0 * (xAB * xAB + yAB * yAB + zAB * zAB) * aA * aB * aP_inv);

	for(int pb = 0; pb < pp; ++pb) {
	  point C = *(ppoints + pb);
	      
	  double xC = (xP - C.x);
	  double yC = (yP - C.y);
	  double zC = (zP - C.z);

	  xPC[pb] = xC;
	  yPC[pb] = yC;
	  zPC[pb] = zC;
	      
	  tval[pb] = aP * (xC * xC + yC * yC + zC * zC);
	}

	for(int pb = 0; pb < pp * nr_roots; ++pb) {
	  *(rts + pb) = 0.0;
	  *(wgh + pb) = 2 * PI * aP_inv * eval * cA * cB;
	}
  
	rys_rw(pp, nr_roots, tval, rts, wgh);  

	int lX = (lB < lA) ? lB : lA;
	int lY = (lB < lA) ? lA : lB;
	int llX = (lB < lA) ? 3 : 3 * (lB + 1) * nr_roots;
	int llY = (lB < lA) ? 3 * (lB + 1) * nr_roots : 3;

	if((lA == 0) && (lB == 0)) {
	  for(int pb = 0; pb < pp; ++pb) {
	    double *int_array = (matrix + shpair_sz * (p + pb));
	    compute_00(beta, int_array, (wgh + 1 * pb));
	  }
	} else if(((lA == 1) && (lB == 0)) || ((lA == 0) && (lB == 1))) {
	  for(int pb = 0; pb < pp; ++pb) {
	    double xC = xPC[pb];
	    double yC = yPC[pb];
	    double zC = zPC[pb];

	    double *int_array = (matrix + shpair_sz * (p + pb));
	    compute_10_01(xPX, yPX, zPX, xC, yC, zC, beta, int_array, (rts + 2 * pb), (wgh + 2 * pb));
	  }
	} else if((lA == 1) && (lB == 1)) {
	  for(int pb = 0; pb < pp; ++pb) {
	    double xC = xPC[pb];
	    double yC = yPC[pb];
	    double zC = zPC[pb];

	    double *int_array = (matrix + shpair_sz * (p + pb));
	    compute_11(xAB, yAB, zAB, xPX, yPX, zPX, xC, yC, zC, aP_inv, beta, int_array, (rts + 2 * pb), (wgh + 2 * pb));
	  }
	} else if(((lA == 2) && (lB == 0)) || ((lA == 0) && (lB == 2))) {
	  for(int pb = 0; pb < pp; ++pb) {
	    double xC = xPC[pb];
	    double yC = yPC[pb];
	    double zC = zPC[pb];

	    double *int_array = (matrix + shpair_sz * (p + pb));
	    compute_20_02(xPX, yPX, zPX, xC, yC, zC, aP_inv, beta, int_array, (rts + 2 * pb), (wgh + 2 * pb));
	  }
	} else {
	  for(int pb = 0; pb < pp; ++pb) {
	    double xC = xPC[pb];
	    double yC = yPC[pb];
	    double zC = zPC[pb];

	    double *int_array = (matrix + shpair_sz * (p + pb));
	    compute_vrr3(nr_roots, 3 * (lB + 1), lX, llX, lY, llY, xPX, yPX, zPX, xC, yC, zC, aP_inv, (rts + nr_roots * pb), (vrr_array + 0), (hrr_array + 0));	    
	    compute_hrr3(nr_roots, 3 * (lB + 1), lX, llX, lY, llY, xAB, yAB, zAB, (vrr_array + 0), (hrr_array + 0));
	    compute_reduction(nr_roots, lA, lB, (wgh + nr_roots * pb), (hrr_array + 0), int_array, beta);
	  }
	}
	    
	beta = 1.0;
      }
    }
  }
  
  free(rts);
  free(wgh);

  free(vrr_array);
  free(hrr_array);
}

#if 0
void compute_integral_shell_pair_pre( int npts,
				      shell_pair shpair, 
				      point *points,
				      double *matrix ) {
  double *rts = (double*) malloc(PB * R_MAX * sizeof(double));
  double *wgh = (double*) malloc(PB * R_MAX * sizeof(double));

  double *vrr_array = (double*) malloc(3 * (Lx + Ly + 1) * R_MAX * sizeof(double));
  double *hrr_array = (double*) malloc(3 * (Lx + 1) * (Ly + 1) * R_MAX * sizeof(double));

  int lA = shpair.lA;
  int lB = shpair.lB;
    
  // nr of roots
  int nr_roots = ((int) ceil((lA + lB) / 2.0)) + 1;

  int value = (lB < lA) ? 1 : -1;    
  double xAB = value * shpair.rAB.x;
  double yAB = value * shpair.rAB.y;
  double zAB = value * shpair.rAB.z;
  
  const int shpair_sz =  (lA+1)*(lA+2) * (lB+1)*(lB+2) / 4;
  for(int p = 0; p < npts; p += PB) {
    int pp = MIN(npts - p, PB);
    point *ppoints = (points + p);
	
    double beta = 0.0;
    prim_pair *prim_pairs = shpair.prim_pairs;
    for(int ij = 0; ij < shpair.nprim_pair; ++ij) { 
      const double aP = prim_pairs[ij].gamma;
      const double aP_inv = 1.0 / aP;

      const double xP = prim_pairs[ij].P.x;
      const double yP = prim_pairs[ij].P.y;
      const double zP = prim_pairs[ij].P.z;
          
      const double xPX = (lB < lA) ? prim_pairs[ij].PA.x : prim_pairs[ij].PB.x;
      const double yPX = (lB < lA) ? prim_pairs[ij].PA.y : prim_pairs[ij].PB.y;
      const double zPX = (lB < lA) ? prim_pairs[ij].PA.z : prim_pairs[ij].PB.z;

      double tval[PB];
      double xPC[PB];
      double yPC[PB];
      double zPC[PB];
	    
      for(int pb = 0; pb < pp; ++pb) {
	point C = *(ppoints + pb);
	      
	double xC = (xP - C.x);
	double yC = (yP - C.y);
	double zC = (zP - C.z);

	xPC[pb] = xC;
	yPC[pb] = yC;
	zPC[pb] = zC;
	      
	tval[pb] = aP * (xC * xC + yC * yC + zC * zC);
      }

      for(int pb = 0; pb < pp * nr_roots; ++pb) {
	*(rts + pb) = 0.0;
	*(wgh + pb) = 2 * PI * aP_inv * prim_pairs[ij].K * prim_pairs[ij].coeff_prod;
      }
  
      rys_rw(pp, nr_roots, tval, rts, wgh);  

      int lX = (lB < lA) ? lB : lA;
      int lY = (lB < lA) ? lA : lB;
      int llX = (lB < lA) ? 3 : 3 * (lB + 1) * nr_roots;
      int llY = (lB < lA) ? 3 * (lB + 1) * nr_roots : 3;

      if((lA == 0) && (lB == 0)) {
	for(int pb = 0; pb < pp; ++pb) {
	  double *int_array = (matrix + shpair_sz * (p + pb));
	  compute_00(beta, int_array, (wgh + 1 * pb));
	}
      } else if(((lA == 1) && (lB == 0)) || ((lA == 0) && (lB == 1))) {
	for(int pb = 0; pb < pp; ++pb) {
	  double xC = xPC[pb];
	  double yC = yPC[pb];
	  double zC = zPC[pb];

	  double *int_array = (matrix + shpair_sz * (p + pb));
	  compute_10_01(xPX, yPX, zPX, xC, yC, zC, beta, int_array, (rts + 2 * pb), (wgh + 2 * pb));
	}
      } else if((lA == 1) && (lB == 1)) {
	for(int pb = 0; pb < pp; ++pb) {
	  double xC = xPC[pb];
	  double yC = yPC[pb];
	  double zC = zPC[pb];

	  double *int_array = (matrix + shpair_sz * (p + pb));
	  compute_11(xAB, yAB, zAB, xPX, yPX, zPX, xC, yC, zC, aP_inv, beta, int_array, (rts + 2 * pb), (wgh + 2 * pb));
	}
      } else if(((lA == 2) && (lB == 0)) || ((lA == 0) && (lB == 2))) {
	for(int pb = 0; pb < pp; ++pb) {
	  double xC = xPC[pb];
	  double yC = yPC[pb];
	  double zC = zPC[pb];

	  double *int_array = (matrix + shpair_sz * (p + pb));
	  compute_20_02(xPX, yPX, zPX, xC, yC, zC, aP_inv, beta, int_array, (rts + 2 * pb), (wgh + 2 * pb));
	}
      } else {
	for(int pb = 0; pb < pp; ++pb) {
	  double xC = xPC[pb];
	  double yC = yPC[pb];
	  double zC = zPC[pb];

	  double *int_array = (matrix + shpair_sz * (p + pb));
	  compute_vrr3(nr_roots, 3 * (lB + 1), lX, llX, lY, llY, xPX, yPX, zPX, xC, yC, zC, aP_inv, (rts + nr_roots * pb), (vrr_array + 0), (hrr_array + 0));	    
	  compute_hrr3(nr_roots, 3 * (lB + 1), lX, llX, lY, llY, xAB, yAB, zAB, (vrr_array + 0), (hrr_array + 0));
	  compute_reduction(nr_roots, lA, lB, (wgh + nr_roots * pb), (hrr_array + 0), int_array, beta);
	}
      }
	    
      beta = 1.0;
    }
  }

  free(rts);
  free(wgh);

  free(vrr_array);
  free(hrr_array);
}
#endif










