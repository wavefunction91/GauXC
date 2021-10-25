#ifndef _MY_BOYS_COMPUTATION
#define _MY_BOYS_COMPUTATION

double boys_asymp(int m, double T);
double boys_reference(int m, double T);
double boys_function(int m, double T);
void boys_function(int m, int npts, const double* T, double* FmT);
void boys_asymp(int npts, int m, const double* T, double* FmT);

#endif
