#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#ifdef _MSC_VER
#include <malloc.h>
#endif

#include "boys.h"

#include "rys_1rw.h"
#include "rys_2rw.h"
#include "rys_3rw.h"
#include "rys_4rw.h"
#include "rys_5rw.h"
#include "rys_xrw.h"

void rys_rw(int nt,
	    int ngqp,
	    double *__restrict tval,
	    double *__restrict rts,
	    double *__restrict wts) {
  switch (ngqp) {
  case 1:
    rys_1rw(nt, tval, rts, wts);
    return;
  case 2:
    rys_2rw(nt, tval, rts, wts);
    return;
  case 3:
    rys_3rw(nt, tval, rts, wts);
    return;
  case 4:
    rys_4rw(nt, tval, rts, wts);
    return;
  case 5:
    rys_5rw(nt, tval, rts, wts);
    return;
  default:
    {
#ifdef _MSC_VER
      double *ryszero = (double *)_malloca(nt * sizeof(double));
#else
      double ryszero[nt];
#endif
      
      for (int n = 0; n < nt; n++) {
	const double t = tval[n];
	if (t == 0.0) {
	  ryszero[n] = 1.0;
	} else if (t <= tmax) {
	  const int tgrid = lround(t * tvstep);
	  const double delta = tgrid * tstep - t;
	  
	  ryszero[n] = (((((boys_table[tgrid][6] * delta * 0.166666666666667 +
			    boys_table[tgrid][5]) * delta * 0.2 +
			   boys_table[tgrid][4]) * delta * 0.25 +
			  boys_table[tgrid][3]) * delta * 0.333333333333333 +
			 boys_table[tgrid][2]) * delta * 0.5 +
			boys_table[tgrid][1]) * delta + boys_table[tgrid][0];
	} else {
	  ryszero[n] = sqrt (3.141592653589793 / t) * .5;
	}
      }
      
      int ntgqp = nt * ngqp;
      int nmom = (ngqp << 1) - 1;
      
      rys_xrw(nt, ntgqp, ngqp, nmom, tval, ryszero, rts, wts);

#ifdef _MSC_VER
      _freea(ryszero);
#endif
      return;
    }
  }
}
