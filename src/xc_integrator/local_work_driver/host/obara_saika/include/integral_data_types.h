#ifndef _MY_INTEGRAL_DATA_TYPES
#define _MY_INTEGRAL_DATA_TYPES
#include <cmath>

typedef struct {
  double x, y, z;
} point;

typedef struct {
  double alpha, coeff;
} coefficients;

typedef struct {
  point origin;
  coefficients *coeff;
  int m, L;
} shells;

typedef struct {
  point P;
  point PA;
  point PB;

  double K;
  double gamma;
  double coeff_prod;
} prim_pair;

typedef struct {
  int lA;
  int lB;
  int nprim_pair;
  point rA, rB;
  point rAB;
  prim_pair* prim_pairs;
} shell_pair;

inline void generate_shell_pair( const shells& A, const shells& B, shell_pair& AB) {
  // L Values
  AB.lA = A.L;
  AB.lB = B.L;

  AB.rA = A.origin;
  AB.rB = B.origin;

  const auto xA = A.origin.x;
  const auto yA = A.origin.y;
  const auto zA = A.origin.z;

  const auto xB = B.origin.x;
  const auto yB = B.origin.y;
  const auto zB = B.origin.z;

  AB.rAB.x = xA - xB;
  AB.rAB.y = yA - yB;
  AB.rAB.z = zA - zB;

  const double dAB = AB.rAB.x*AB.rAB.x + AB.rAB.y*AB.rAB.y + AB.rAB.z*AB.rAB.z;
  
  const int nprim_A = A.m;
  const int nprim_B = B.m;
  const int np = nprim_A * nprim_B;

  AB.nprim_pair = np;
  AB.prim_pairs = new prim_pair[np];
  for( int i = 0, ij = 0; i < nprim_A; ++i       )
  for( int j = 0        ; j < nprim_B; ++j, ++ij ) {
    auto& pair = AB.prim_pairs[ij];
    pair.coeff_prod = A.coeff[i].coeff * B.coeff[j].coeff;

    const auto alpha_A = A.coeff[i].alpha;
    const auto alpha_B = B.coeff[j].alpha;

    pair.gamma = alpha_A + alpha_B;
    const auto gamma_inv = 1. / pair.gamma;

    pair.P.x = (alpha_A * xA + alpha_B * xB) * gamma_inv;
    pair.P.y = (alpha_A * yA + alpha_B * yB) * gamma_inv;
    pair.P.z = (alpha_A * zA + alpha_B * zB) * gamma_inv;

    pair.PA.x = pair.P.x - xA;
    pair.PA.y = pair.P.y - yA;
    pair.PA.z = pair.P.z - zA;

    pair.PB.x = pair.P.x - xB;
    pair.PB.y = pair.P.y - yB;
    pair.PB.z = pair.P.z - zB;

    pair.K = 2 * M_PI * gamma_inv *
      std::exp( - alpha_A * alpha_B * dAB * gamma_inv );
  }
}

#endif
