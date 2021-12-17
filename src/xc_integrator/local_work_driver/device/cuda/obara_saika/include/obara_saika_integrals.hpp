#ifndef __MY_INTEGRAL_OBARA_SAIKA
#define __MY_INTEGRAL_OBARA_SAIKA

void compute_integral_shell_pair(size_t npts,
				 int i,
				 int j,
				 shells *shell_list,
				 double *points,
				 double *Xi,
				 double *Xj,
				 int ldX,
				 double *Gi,
				 double *Gj,
				 int ldG, 
				 double *weights,
				 double *boys_table);

#endif
