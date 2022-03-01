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
  point rA, rB;
  point rAB;
} shell_pair;

#endif