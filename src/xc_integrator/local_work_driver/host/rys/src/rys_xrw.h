#ifndef RYS_XRW_H_
#define RYS_XRW_H_

void rys_xrw(int nt,
	     int ntgqp,
	     int ngqp,
	     int nmom,
	     const double tval[restrict],
	     const double ryszero[restrict],
	     double rts[restrict],
	     double wts[restrict]);

#endif
