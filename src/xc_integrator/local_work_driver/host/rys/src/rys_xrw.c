#include <stdint.h>
#include <stddef.h>
#include <math.h>
#include <assert.h>
#include "jacobi.h"

void rys_xrw(int nt,
	      int ntgqp,
	      int ngqp,
	      int nmom,
	      const double tval[restrict],
	      const double ryszero[restrict],
	      double rts[restrict],
	      double wts[restrict]) {
  double a[nmom];
  double b[nmom-1];
  double mom[nmom];
  double dia[ngqp];
  double off[ngqp];
  double row1[nmom];
  double row2[nmom];

  int nrts = 0;
  for (int n = 0; n < nt; n += 1) {
    const double t = tval[n];
    const double momzero = ryszero[n];
    if (t <= 15.0) {
      
      assert(nmom <= 30);

      if (t <= 1.0e-16) {
	const int imax = (nmom < 16) ? nmom : 16;
	a[0] = ajac[0];
	mom[0] = csmall[0] * t;
	double tpower = t;
	for (int i = 2; i <= imax; ++i) {
	  tpower *= t;
	  a[i-1] = ajac[i - 1];
	  b[i-2] = bjac[i - 2];
	  mom[i-1] = csmall[i - 1] * tpower;
	}
	for (int i = imax + 1; i <= nmom; ++i) {
	  a[i-1] = ajac[i - 1];
	  b[i-2] = bjac[i - 2];
	  mom[i-1] = 0.;
	}
      } else {
	int imax;
	if (nmom <= 5) {
	  if (t < 1.0e-6) {
	    imax = nmom + 1;
	  } else if (t < .1) {
	    imax = nmom + 3;
	  } else if (t < 2.) {
	    imax = nmom + 7;
	  } else if (t < 10.) {
	    imax = nmom + 13;
	  } else {
	    imax = nmom + 22;
	  }
	} else {
	  if (t < 1.0e-6) {
	    imax = nmom;
	  } else if (t < .1) {
	    imax = nmom + 2;
	  } else if (t < 2.) {
	    imax = nmom + 4;
	  } else if (t < 10.) {
	    imax = nmom + 8;
	  } else {
	    imax = nmom + 16;
	  }
	}

	double momi = 1.0e-300;
	double momip1 = 0.0;
	const double tinvhf = .5 / t;
	double r1 = (double) ((imax << 1) + 5);
	for (int i = imax + 1; i >= nmom + 2; --i) {
	  r1 -= 2.0;
	  const double r = r1 * tinvhf + r2[i - 1];
	  const double momim1 = sinv[i - 1] * (momip1 - r * momi);
	  momip1 = momi;
	  momi = momim1;
	}
	for (int i = nmom + 1; i >= 2; --i) {
	  r1 -= 2.0;
	  const double r = r1 * tinvhf + r2[i - 1];
	  const double momim1 = sinv[i - 1] * (momip1 - r * momi);
	  mom[i - 2] = momim1;
	  momip1 = momi;
	  momi = momim1;
	}

	const double r = tinvhf * 3.0 + r2[0];
	const double zmom = sinv[0] * (momip1 - r * momi);
	assert(fabs(zmom) >= 1.0e-300);
	a[0] = ajac[0];
	const double zinv = 1. / zmom;
	mom[0] *= zinv;
	for (int i = 2; i <= nmom; ++i) {
	  a[i-1] = ajac[i - 1];
	  b[i-2] = bjac[i - 2];
	  mom[i-1] *= zinv;
	}
      }
    } else {
      const double texp = exp(-t);
      const double tinv = 1.0 / t;
      const double tinv2 = tinv * 2.;
      const double tinvhf = tinv * .5;
      const double tinvsq = tinv * tinv;
      const double scale = -tinvhf * texp / momzero;
      if (nmom == 1) {
	a[0] = tinvhf;
	mom[0] = scale;
      } else {
	a[0] = tinvhf;
	a[1] = tinvhf + tinv2;
	b[0] = tinvsq * .5;
	mom[0] = scale;
	double r = 1. - tinv * 1.5;
	mom[1] = scale * r;
	double s = 0.0;
	double binc = 0.5;
	double sinc = -0.5;
	double lim2 = r;
	double lim3 = 1.0;
	for (int i = 3; i <= nmom; ++i) {
	  binc += 2.;
	  a[i-1] = a[i-2] + tinv2;
	  b[i-2] = b[i - 3] + binc * tinvsq;
	  sinc += 2.;
	  r -= tinv2;
	  s += sinc * tinvsq;
	  const double lim1 = r * lim2 - s * lim3;
	  mom[i-1] = scale * lim1;
	  lim3 = lim2;
	  lim2 = lim1;
	}
      }
    }

    if (ngqp == 1) {
      dia[0] = mom[0] + a[0];
    } else if (ngqp == 2) {
      const double sigma = mom[0] + a[0];
      dia[0] = sigma;
      const double theta = (a[1] - sigma) * mom[0] + mom[1] + b[0];
      off[0] = sqrt(theta);
      dia[1] = ((a[2] - sigma) * mom[1] + mom[2] + b[1] * mom[0]) / theta - mom[0] + a[1];
    } else {
      const int imax = ngqp - 1;
      static int jmax = 0;
      jmax = ngqp + imax;
      for (int j = 1; j <= jmax; ++j) {
	row1[j-1] = mom[j-1];
      }
      double sigma = row1[0] + a[0];
      dia[0] = sigma;


      row2[0] = (a[1] - sigma) * row1[0] + row1[1] + b[0];
      double theta = row2[0];
      off[0] = sqrt(theta);
      --jmax;
      for (int j = 2; j <= jmax; ++j) {
	row2[j-1] = (a[j] - sigma) * row1[j-1] + row1[j] + b[j-1] * row1[j - 2];
      }
      sigma = row2[1] / theta - row1[0] + a[1];
      dia[1] = sigma;

      for (int i = 2; i <= imax; ++i) {
	--jmax;
	if (i % 2 == 0) {
	  for (int j = i; j <= jmax; ++j) {
	    row1[j-1] = (a[j] - sigma) * row2[j-1] + row2[j] + b[j-1] * row2[j - 2] - theta * row1[j-1];
	  }
	  sigma = a[i] - row2[i-1] / row2[i - 2] + row1[i] / row1[i-1];
	  theta = row1[i-1] / row2[i - 2];
	} else {
	  for (int j = i; j <= jmax; ++j) {
	    row2[j-1] = (a[j] - sigma) * row1[j-1] + row1[j] + b[j-1] * row1[j - 2] - theta * row2[j-1];
	  }
	  sigma = a[i] - row1[i-1] / row1[i - 2] + row2[i] / row2[i-1];
	  theta = row2[i-1] / row1[i - 2];
	}
	dia[i] = sigma;
	off[i-1] = sqrt(theta);
      }
    }

    if (ngqp == 1) {
      ++nrts;
      rts[nrts-1] = dia[0];
      wts[nrts-1] *= momzero;
    } else {
      a[0] = 1.0;
      for (int j = 2; j <= ngqp; ++j) {
	a[j-1] = 0.0;
      }

      off[ngqp-1] = 0.0;
      int m, iter = 0;
      for (int j = 1; j <= ngqp; ++j) {
      next_iteration:
	for (m = j; m < ngqp; ++m) {
	  const double test1 = fabs(dia[m-1]) + fabs(dia[m]);
	  const double test2 = test1 + fabs(off[m-1]);
	  if (test2 == test1) {
	    break;
	  }
	}
	double p = dia[j-1];
	if (m != j) {
	  assert(iter != 30);
	  ++iter;
	  double g = (dia[j] - p) / (off[j-1] * 2.);
	  double r = sqrt(g * g + 1.);
	  g = dia[m-1] - p + off[j-1] / (g + copysign(r, g));
	  double s = 1.0;
	  double c = 1.0;
	  p = 0.0;
	  for (int i = m - 1; i >= j; --i) {
	    double f = s * off[i-1];
	    const double d = c * off[i-1];
	    r = sqrt(f * f + g * g);
	    off[i] = r;
	    if (r == 0.0) {
	      dia[i] -= p;
	      off[m-1] = 0.0;
	      goto next_iteration;
	    }
	    s = f / r;
	    c = g / r;
	    g = dia[i] - p;
	    r = (dia[i-1] - g) * s + c * 2. * d;
	    p = s * r;
	    dia[i] = g + p;
	    g = c * r - d;
	    f = a[i];
	    a[i] = s * a[i-1] + c * f;
	    a[i-1] = c * a[i-1] - s * f;
	  }
	  dia[j-1] -= p;
	  off[j-1] = g;
	  off[m-1] = 0.0;
	  goto next_iteration;
	}
      }

      for (int i = 1; i <= ngqp; ++i) {
	const double root = dia[i-1];

	assert((root >= 0.0) && (root <= 1.0));
	rts[nrts+i-1] = root;

	const double ai = a[i-1];
	wts[nrts+i-1] *= momzero * (ai * ai);
      }
      nrts += ngqp;
    }
  }
}
