#include <gauxc/xc_integrator/default_xc_host_integrator.hpp>

#include "host_weights.hpp"
#include "host_collocation.hpp"
#include "integrator_common.hpp"
#include "blas.hpp"

namespace GauXC  {
namespace detail {


template <typename _F1, typename _F2>
void submat_set(int32_t M, int32_t N, int32_t MSub, 
  int32_t NSub, _F1 *ABig, int32_t LDAB, _F2 *ASmall, 
  int32_t LDAS, 
  std::vector<std::pair<int32_t,int32_t>> &submat_map) {


  int32_t i(0);
  for( auto& iCut : submat_map ) {
    int32_t deltaI = iCut.second - iCut.first;
    int32_t j(0);
  for( auto& jCut : submat_map ) {
    int32_t deltaJ = jCut.second - jCut.first;
  
    auto* ABig_use   = ABig   + iCut.first + jCut.first * LDAB;
    auto* ASmall_use = ASmall + i          + j          * LDAS;


    GauXC::blas::lacpy( 'A', deltaI, deltaJ, ABig_use, LDAB, 
                         ASmall_use, LDAS );

  
    j += deltaJ;
  }
    i += deltaI;
  }
  

}

template <typename F, size_t n_deriv>
void process_batches_host_replicated_p(
  XCWeightAlg            weight_alg,
  const functional_type& func,
  const BasisSet<F>&     basis,
  const Molecule   &     mol,
  const MolMeta    &     meta,
  XCHostData<F>    &     host_data,
  std::vector< XCTask >& tasks,
  const F*               P,
  F*                     VXC,
  F*                     exc,
  F*                     n_el
) {

  partition_weights_host( weight_alg, mol, meta, tasks );

  size_t ntasks = tasks.size();
  for( size_t iT = 0; iT < ntasks; ++iT ) {

    auto& task = tasks[iT];

    const size_t  npts    = task.points.size();
    const size_t  nbe     = task.nbe;
    const size_t  nshells = task.shell_list.size();

    const double* points      = task.points.data()->data();
    const double* weights     = task.weights.data();
    const int32_t* shell_list = task.shell_list.data();

    double* basis_eval = host_data.basis_eval.data();
    double* den_eval   = host_data.den_scr.data();
    double* nbe_scr    = host_data.nbe_scr.data();
    double* zmat       = host_data.zmat.data();

    double* eps        = host_data.eps.data();
    double* gamma      = host_data.gamma.data();
    double* vrho       = host_data.vrho.data();
    double* vgamma     = host_data.vgamma.data();

    double* dbasis_x_eval = nullptr;
    double* dbasis_y_eval = nullptr;
    double* dbasis_z_eval = nullptr;
    double* dden_x_eval = nullptr;
    double* dden_y_eval = nullptr;
    double* dden_z_eval = nullptr;

    if( n_deriv > 0 ) {
      dbasis_x_eval = basis_eval    + npts * nbe;
      dbasis_y_eval = dbasis_x_eval + npts * nbe;
      dbasis_z_eval = dbasis_y_eval + npts * nbe;
      dden_x_eval   = den_eval    + npts;
      dden_y_eval   = dden_x_eval + npts;
      dden_z_eval   = dden_y_eval + npts;
    }


    // Get the submatrix map for batch
    auto submat_map = gen_compressed_submat_map( basis, task.shell_list );


    // Evaluate Collocation Matrix 
    if( n_deriv == 1 )
      eval_collocation_deriv1( npts, nshells, nbe, points, basis, shell_list, 
                               basis_eval, dbasis_x_eval, dbasis_y_eval, 
                               dbasis_z_eval );
    else
      eval_collocation( npts, nshells, nbe, points, basis, shell_list, basis_eval );


    // Extrat Submatrix
    const F* den_ptr_use = P;
    if( nbe != basis.nbf() ) {
      submat_set( basis.nbf(), basis.nbf(), nbe, nbe, P, basis.nbf(),
                  nbe_scr, nbe, submat_map );
      den_ptr_use = nbe_scr;
    } 

    // Z = P * BF
    GauXC::blas::gemm( 'N', 'N', nbe, npts, nbe, 1., den_ptr_use, nbe,
                       basis_eval, nbe, 0., zmat, nbe );
    

    // Evaluate the density 
    for( int32_t i = 0; i < npts; ++i ) {

      const size_t ioff = size_t(i) * nbe;
      const F*     zmat_i = zmat + ioff;

      den_eval[i] = 
        2. * GauXC::blas::dot( nbe, basis_eval + ioff, 1, zmat_i, 1 );

      if( n_deriv > 0 ) {
        const double dx = 
          4. * GauXC::blas::dot( nbe, dbasis_x_eval + ioff, 1, zmat_i, 1 );
        const double dy = 
          4. * GauXC::blas::dot( nbe, dbasis_y_eval + ioff, 1, zmat_i, 1 );
        const double dz = 
          4. * GauXC::blas::dot( nbe, dbasis_z_eval + ioff, 1, zmat_i, 1 );

        dden_x_eval[i] = dx;
        dden_y_eval[i] = dy;
        dden_z_eval[i] = dz;

        gamma[i] = dx*dx + dy*dy + dz*dz;
      }

    }


    // Evaluate XC functional
    if( func.is_gga() )
      func.eval_exc_vxc( npts, den_eval, gamma, eps, vrho, vgamma );
    else
      func.eval_exc_vxc( npts, den_eval, eps, vrho );

    // Scalar integrations
    if( n_el )
      for( int32_t i = 0; i < npts; ++i ) *n_el += weights[i] * den_eval[i];

    for( int32_t i = 0; i < npts; ++i ) *exc += weights[i] * den_eval[i] * eps[i];
    
  }

  //std::cout << std::scientific << std::setprecision(12);
  //std::cout << "NEL = " << *n_el << std::endl;
}


#define HOST_IMPL( F, ND ) \
template \
void process_batches_host_replicated_p<F, ND>(\
  XCWeightAlg            weight_alg,\
  const functional_type& func,\
  const BasisSet<F>&     basis,\
  const Molecule   &     mol,\
  const MolMeta    &     meta,\
  XCHostData<F>    &     host_data,\
  std::vector< XCTask >& local_work,\
  const F*               P,\
  F*                     VXC,\
  F*                     exc,\
  F*                     n_el\
) 

HOST_IMPL( double, 0 );
HOST_IMPL( double, 1 );

}
}
