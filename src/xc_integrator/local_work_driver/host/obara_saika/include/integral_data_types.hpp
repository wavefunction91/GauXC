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

  double K_coeff_prod;
  double gamma;
} prim_pair;

typedef struct {
  int lA;
  int lB;
  int nprim_pair;
  point rA, rB;
  point rAB;
  prim_pair* prim_pairs;
} shell_pair;

#endif
