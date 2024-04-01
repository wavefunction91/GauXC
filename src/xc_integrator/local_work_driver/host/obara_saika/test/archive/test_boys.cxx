/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <cmath>

#include <sys/time.h>

#include "../include/chebyshev_boys_computation.hpp"

int main(int argc, char **argv) {

  int runs = atoi(argv[1]);
  double TVAL = atof(argv[2]);
  
  GauXC::gauxc_boys_init();

  struct timeval t0, t1;
  
  double f0, f1, f2, f3, f4, f5, f6, f7, f8;
  double g0, g1, g2, g3, g4, g5, g6, g7, g8;
  
  long long sum0 = 0;
  for(int r = 0; r < runs; ++r) {
    gettimeofday(&t0, NULL);
    f0 = GauXC::gauxc_boys_element<0>(TVAL);
    f1 = GauXC::gauxc_boys_element<1>(TVAL);
    f2 = GauXC::gauxc_boys_element<2>(TVAL);
    f3 = GauXC::gauxc_boys_element<3>(TVAL);
    f4 = GauXC::gauxc_boys_element<4>(TVAL);
    f5 = GauXC::gauxc_boys_element<5>(TVAL);
    f6 = GauXC::gauxc_boys_element<6>(TVAL);
    f7 = GauXC::gauxc_boys_element<7>(TVAL);
    f8 = GauXC::gauxc_boys_element<8>(TVAL);
    gettimeofday(&t1, NULL);

    sum0 += (t1.tv_sec-t0.tv_sec)*1000000LL + t1.tv_usec-t0.tv_usec;
  }

  long long sum1 = 0;
  for(int r = 0; r < runs; ++r) {
    gettimeofday(&t0, NULL);
    double e_TVAL_neg = std::exp(-TVAL);
    double TVAL_rec = 1 / (2 * TVAL);

    e_TVAL_neg = (-1.0) * e_TVAL_neg * TVAL_rec;
  
    g0 = GauXC::gauxc_boys_element<0>(TVAL);

    g1 =  1.0 * g0 * TVAL_rec + e_TVAL_neg;
    g2 =  3.0 * g1 * TVAL_rec + e_TVAL_neg;
    g3 =  5.0 * g2 * TVAL_rec + e_TVAL_neg;
    g4 =  7.0 * g3 * TVAL_rec + e_TVAL_neg;
    g5 =  9.0 * g4 * TVAL_rec + e_TVAL_neg;
    g6 = 11.0 * g5 * TVAL_rec + e_TVAL_neg;
    g7 = 13.0 * g6 * TVAL_rec + e_TVAL_neg;
    g8 = 15.0 * g7 * TVAL_rec + e_TVAL_neg;
    gettimeofday(&t1, NULL);
    
    sum1 += (t1.tv_sec-t0.tv_sec)*1000000LL + t1.tv_usec-t0.tv_usec;
  }

  printf("%lf - %lf = %e\n", f0, g0, f0 - g0);
  printf("%lf - %lf = %e\n", f1, g1, f1 - g1);
  printf("%lf - %lf = %e\n", f2, g2, f2 - g2);
  printf("%lf - %lf = %e\n", f3, g3, f3 - g3);
  printf("%lf - %lf = %e\n", f4, g4, f4 - g4);
  printf("%lf - %lf = %e\n", f5, g5, f5 - g5);
  printf("%lf - %lf = %e\n", f6, g6, f6 - g6);
  printf("%lf - %lf = %e\n", f7, g7, f7 - g7);
  printf("%lf - %lf = %e\n", f8, g8, f8 - g8);

  printf("%lf\t%lf\n", sum0 / ((double) (1.0 * runs)), sum1 / ((double) (1.0 * runs)));
  
  GauXC::gauxc_boys_finalize();
  
  return 0;
}
