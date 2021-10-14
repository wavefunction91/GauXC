#ifndef _MY_INTEGRAL_DATA_TYPES
#define _MY_INTEGRAL_DATA_TYPES

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

#endif
