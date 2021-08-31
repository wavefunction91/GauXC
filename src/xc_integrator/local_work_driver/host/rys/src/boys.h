#ifndef BOYS_H_
#define BOYS_H_

#include <stddef.h>

#define NGRID    920
#define MGRID     10

static const double tmax = 46.0;
static const double tvstep = 20.0;
static const double tstep = 0.05;

extern double boys_table[NGRID + 1][MGRID + 1];

#endif
