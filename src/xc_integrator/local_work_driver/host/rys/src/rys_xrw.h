#ifndef RYS_XRW_H_
#define RYS_XRW_H_

void rys_xrw(int nt,
	     int ntgqp,
	     int ngqp,
	     int nmom,
	     const double *__restrict tval,
	     const double *__restrict ryszero,
	     double *__restrict rts,
	     double *__restrict wts);

#endif
