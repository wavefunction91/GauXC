#ifndef __MY_INTEGRAL_OBARA_SAIKA
#define __MY_INTEGRAL_OBARA_SAIKA

#include "integral_data_types.hpp"

namespace XCPU {
<<<<<<< HEAD

  void generate_shell_pair( const shells& A, const shells& B, shell_pair& AB);
  
  void compute_integral_shell_pair(size_t npts,
				   int i,
				   int j,
				   shells *shell_list,
=======
  void generate_shell_pair( const shells& A, const shells& B, shell_pair& AB);
  void compute_integral_shell_pair(size_t npts,
				   int is_diag,
				   int lA,
				   int lB,
				   shell_pair &shpair,
>>>>>>> 95003a054eaeb7ccfeabc7e81aaa0439f60d6524
				   double *points,
				   double *Xi,
				   double *Xj,
				   int ldX,
				   double *Gi,
				   double *Gj,
				   int ldG, 
				   double *weights, 
				   double *boys_table);
<<<<<<< HEAD

  void compute_integral_shell_pair_v0(size_t npts,
				      int is_diag,
				      int lA,
				      int lB,
				      shell_pair &shpair,
				      double *points,
				      double *Xi,
				      double *Xj,
				      int ldX,
				      double *Gi,
				      double *Gj,
				      int ldG, 
				      double *weights, 
				      double *boys_table);
=======
>>>>>>> 95003a054eaeb7ccfeabc7e81aaa0439f60d6524
}

#endif
