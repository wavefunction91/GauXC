#ifndef JACOBI_H_
#define JACOBI_H_

/* (2i+1) / ((4i-1) * (4i+3)) for i = 1...100 */
extern const double r2[100];

/* ((4i-3) * (4i+1) * square(4i-1)) / (2i * (2i+1) * square(2i-1)) for i = 1...100 */
extern const double sinv[100];

extern const double csmall[16];

/* (4i * (2i+1) - 1) / ((4i+3) * (4i-1)) for i = 0...99 */
extern const double ajac[100];

/* (4*square(i) * (4*square(i) - 4*i + 1)) / ((4i-3) * (4i+1) * square(4i-1)) for i = 1...99 */
extern const double bjac[99];

#endif // JACOBI_H_
