#ifndef _MY_SIMD_INSTRUCTIONS
#define _MY_SIMD_INSTRUCTIONS

#define NPTS_LOCAL 64

#define SCALAR_TYPE double

#define SCALAR_LENGTH 1

#define SCALAR_SET1(x) x
#define SCALAR_ZERO() 0.0

#define SCALAR_LOAD(x) *(x)
#define SCALAR_STORE(x, y) *(x) = y

#define SCALAR_ADD(x, y) x + y
#define SCALAR_SUB(x, y) x - y

#define SCALAR_MUL(x, y) x * y
#define SCALAR_FMA(x, y, z) z + x * y 
#define SCALAR_FNMA(x, y, z) z - x * y

#define SCALAR_RECIPROCAL(x) 1.0 / (1.0 * x)

#endif
