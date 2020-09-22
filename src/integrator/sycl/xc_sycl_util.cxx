#include <gauxc/xc_integrator/xc_sycl_util.hpp>
#include <gauxc/util/sycl_util.hpp>

#include "sycl_weights.hpp"
#include "collocation_device.hpp"
#include "sycl_pack_density.hpp"
#include "sycl_inc_potential.hpp"
#include "sycl_eval_denvars.hpp"
#include "sycl_zmat.hpp"
#include "integrator_common.hpp"

#include "mklsycl_extensions.hpp"

#include <Eigen/Core>

inline static void unused() { }
template <typename T, typename... Args>
inline static void unused( const T& t, Args&&... args ) {
  (void)(t);
  unused( std::forward<Args>(args)... );
}

namespace GauXC  {
namespace integrator::sycl {

using namespace GauXC::sycl::blas;

using host_task_iterator = std::vector<XCTask>::iterator;

template <typename F>
using sycl_task_iterator = typename std::vector<XCTaskDevice<F>>::iterator;

template <typename F, size_t n_deriv>
void process_batches_sycl_replicated_all_device(
  XCWeightAlg            weight_alg,
  const functional_type& func,
  XCSyclData<F>&         sycl_data,
  sycl_task_iterator<F>  task_begin,
  sycl_task_iterator<F>  task_end) {

  const auto ntasks = std::distance( task_begin, task_end );
  const auto nbf    = sycl_data.nbf;

  // Get batch statistics for batches to process
  auto nbe_comparator =
    []( const auto& a, const auto& b ){ return a.nbe < b.nbe; };
  auto npts_comparator =
    []( const auto& a, const auto& b ){ return a.npts < b.npts; };
  auto nshells_comparator =
    []( const auto& a, const auto& b ){ return a.nshells < b.nshells; };

  auto [min_nbe_it, max_nbe_it] =
    std::minmax_element( task_begin, task_end, nbe_comparator );
  auto [min_npts_it, max_npts_it] =
    std::minmax_element( task_begin, task_end, npts_comparator );
  auto [min_nshells_it, max_nshells_it] =
    std::minmax_element( task_begin, task_end, nshells_comparator );

  const auto min_nbe     = min_nbe_it->nbe;
  const auto max_nbe     = max_nbe_it->nbe;
  const auto min_npts    = min_npts_it->npts;
  const auto max_npts    = max_npts_it->npts;
  const auto min_nshells = min_nshells_it->nshells;
  const auto max_nshells = max_nshells_it->nshells;

  unused( min_nbe, min_npts, min_nshells );

  const size_t total_npts =
    std::accumulate( task_begin, task_end, 0ul,
                     []( const auto& a, const auto& b ) { return a + b.npts; } );


  // Aliases
  cl::sycl::queue syclQue = *(sycl_data.master_queue);

  auto* dmat_device         = sycl_data.dmat_device;

  auto* shells_device       = sycl_data.shells_device;
  auto* tasks_device        = sycl_data.device_tasks;
  auto* dmat_array_device   = sycl_data.dmat_array_device;
  auto* zmat_array_device   = sycl_data.zmat_array_device;
  auto* bf_array_device     = sycl_data.bf_array_device;
  auto* weights_device      = sycl_data.weights_device_buffer;
  auto* dist_scratch_device = sycl_data.dist_scratch_device;

  auto* den_eval_device     = sycl_data.den_eval_device;
  auto* dden_x_eval_device  = sycl_data.den_x_eval_device;
  auto* dden_y_eval_device  = sycl_data.den_y_eval_device;
  auto* dden_z_eval_device  = sycl_data.den_z_eval_device;

  auto* eps_eval_device     = sycl_data.eps_eval_device;
  auto* gamma_eval_device   = sycl_data.gamma_eval_device;
  auto* vrho_eval_device    = sycl_data.vrho_eval_device;
  auto* vgamma_eval_device  = sycl_data.vgamma_eval_device;


  auto* exc_device     = sycl_data.exc_device;
  auto* vxc_device     = sycl_data.vxc_device;
  auto* nel_device     = sycl_data.nel_device;
  auto* acc_scr_device = sycl_data.acc_scr_device;

  auto* m_array_device      = sycl_data.m_array_device;
  auto* n_array_device      = sycl_data.n_array_device;
  auto* k_array_device      = sycl_data.k_array_device;
  auto* lda_array_device    = sycl_data.lda_array_device;
  auto* ldb_array_device    = sycl_data.ldb_array_device;
  auto* ldc_array_device    = sycl_data.ldc_array_device;

  const auto* rab_device          = sycl_data.rab_device;
  const auto* coords_device       = sycl_data.coords_device;
  const auto* points_device       = sycl_data.points_device_buffer;
  const auto* iparent_device      = sycl_data.iparent_device_buffer;
  const auto* dist_nearest_device = sycl_data.dist_nearest_buffer;


  // Evaluate Partition Weights
  partition_weights_sycl_SoA( weight_alg, total_npts, sycl_data.natoms,
                              points_device, iparent_device, dist_nearest_device,
                              rab_device, coords_device, weights_device,
                              dist_scratch_device, &syclQue );

  std::cout << "AFTER WEIGHTS" << std::endl;

  // Evaluate Collocation
  if constexpr ( n_deriv == 1 )
    eval_collocation_masked_combined_deriv1( ntasks, max_npts, max_nshells,
                                             shells_device, tasks_device,
                                             &syclQue );
  else
    eval_collocation_masked_combined( ntasks, max_npts, max_nshells, shells_device,
                                      tasks_device, &syclQue );

  std::cout << "AFTER COLLOCATION" << std::endl;
  // Pack Density Submatrices
  task_pack_density_matrix( ntasks, tasks_device, dmat_device, nbf, &syclQue );

  std::cout << "AFTER PACK" << std::endl;

  // Form Z = P * X
  for( auto it = task_begin; it != task_end; ++it )
    oneapi::mkl::blas::column_major::gemm( syclQue, oneapi::mkl::transpose::nontrans, oneapi::mkl::transpose::nontrans,
		                           it->nbe, it->npts, it->nbe, 1., it->nbe_scr, it->nbe, it->bf, it->nbe,
					   0., it->zmat, it->nbe );

  std::cout << "AFTER GEMM" << std::endl;

  // Zero UVars
  util::sycl_set_zero_async( total_npts, den_eval_device, syclQue, "DenZero" );
  if( func.is_gga() ) {
    util::sycl_set_zero_async( total_npts, dden_x_eval_device, syclQue,
                               "DenXZero" );
    util::sycl_set_zero_async( total_npts, dden_y_eval_device, syclQue,
                               "DenYZero" );
    util::sycl_set_zero_async( total_npts, dden_z_eval_device, syclQue,
                               "DenZZero" );
  }

  // Evaluate UVars
  if( func.is_gga() ) {
    eval_uvars_gga_device( ntasks, max_nbe, max_npts, tasks_device, &syclQue );
    eval_vvars_gga_device( total_npts, dden_x_eval_device, dden_y_eval_device,
                           dden_z_eval_device, gamma_eval_device, &syclQue );
  } else {
    eval_uvars_lda_device( ntasks, max_nbe, max_npts, tasks_device, &syclQue );
  }

  std::cout << "AFTER UVARS" << std::endl;
  // Evaluate XC Functional
  if( func.is_gga() )
    func.eval_exc_vxc_device( total_npts, den_eval_device, gamma_eval_device,
                              eps_eval_device, vrho_eval_device,
                              vgamma_eval_device, &syclQue );
  else
    func.eval_exc_vxc_device( total_npts, den_eval_device, eps_eval_device,
                              vrho_eval_device, &syclQue );

  std::cout << "AFTER XC" << std::endl;
  // Factor weights into XC output
  hadamard_product( &syclQue, total_npts, 1, weights_device, 1,
                    eps_eval_device, 1 );
  hadamard_product( &syclQue, total_npts, 1, weights_device, 1,
                    vrho_eval_device, 1 );
  if( func.is_gga() )
    hadamard_product( &syclQue, total_npts, 1, weights_device, 1,
                      vgamma_eval_device, 1 );

  // Accumulate EXC / NEL
  gdot( &syclQue, total_npts, weights_device, 1,
        den_eval_device, 1, acc_scr_device, nel_device );
  gdot( &syclQue, total_npts, eps_eval_device, 1,
        den_eval_device, 1, acc_scr_device, exc_device );

  std::cout << "AFTER SCALAR" << std::endl;
  // Evaluate Z Matrix
  if( func.is_gga() )
    zmat_gga_sycl( ntasks, max_nbe, max_npts, tasks_device, &syclQue );
  else
    zmat_lda_sycl( ntasks, max_nbe, max_npts, tasks_device, &syclQue );

  std::cout << "AFTER ZMAT" << std::endl;

  // // Accumulate packed VXC = X * Z**T + Z * X**T
  // // XXX: Only updates LT
  for( auto it = task_begin; it != task_end; ++it )
    oneapi::mkl::blas::column_major::syr2k( syclQue,
                                            oneapi::mkl::uplo::lower, oneapi::mkl::transpose::nontrans,
                                            it->nbe, it->npts,
                                            1., it->bf, it->nbe,
                                            it->zmat, it->nbe,
                                            0., it->nbe_scr, it->nbe );

  std::cout << "AFTER SYR2K" << std::endl;

  // Increment global VXC
  task_inc_potential( ntasks, tasks_device, vxc_device, nbf, &syclQue );


  // Synchronize on master queue
  // XXX: There's no lifetime issues in this driver, should look into
  //      avoid this sync to allow for overlap with the host packing
  util::sycl_device_sync( syclQue );

}


template <typename F, size_t n_deriv>
void process_batches_sycl_replicated_p(
  XCWeightAlg            weight_alg,
  const functional_type& func,
  const BasisSet<F>&     basis,
  const Molecule   &     mol,
  const MolMeta    &     meta,
  XCSyclData<F>    &     sycl_data,
  std::vector< XCTask >& tasks,
  const F*               P,
  F*                     VXC,
  F*                     EXC,
  F*                     NEL) {

  auto task_comparator = []( const XCTask& a, const XCTask& b ) {
    return (a.points.size() * a.nbe) > (b.points.size() * b.nbe);
  };
  std::sort( tasks.begin(), tasks.end(), task_comparator );

  cl::sycl::queue syclQue = *(sycl_data.master_queue);

  const auto nbf    = basis.nbf();
  const auto natoms = meta.natoms();

  // Send static data to the device

  // Density
  if( not sycl_data.denpack_host )
      util::sycl_copy( nbf * nbf, sycl_data.dmat_device, P, syclQue, "P H2D" );

  // Shells: TODO avoid host copy?
  std::vector<Shell<F>> shells( basis );
  util::sycl_copy( shells.size(), sycl_data.shells_device, shells.data(),
                   syclQue, "Shells H2D" );

  // RAB
  util::sycl_copy( natoms * natoms, sycl_data.rab_device, meta.rab().data(),
                   syclQue, "RAB H2D" );

  // Atomic coordinates
  std::vector<double> coords( 3*natoms );
  for( auto i = 0ul; i < natoms; ++i ) {
    coords[ 3*i + 0 ] = mol[i].x;
    coords[ 3*i + 1 ] = mol[i].y;
    coords[ 3*i + 2 ] = mol[i].z;
  }
  util::sycl_copy( 3 * natoms, sycl_data.coords_device, coords.data(),
                   syclQue, "Coords H2D" );


  // Zero out XC quantities
  util::sycl_set_zero( nbf * nbf, sycl_data.vxc_device, syclQue, "VXC Zero" );
  util::sycl_set_zero( 1        , sycl_data.exc_device, syclQue, "EXC Zero" );
  util::sycl_set_zero( 1        , sycl_data.nel_device, syclQue, "NEL Zero" );

  // Processes batches in groups that saturadate available device memory
  auto task_it = tasks.begin();
  while( task_it != tasks.end() ) {

    // Determine next task batch, send relevant data to device
    auto [it, tasks_device] =
      sycl_data.generate_buffers( basis, task_it, tasks.end() );

    // Process the batches
    process_batches_sycl_replicated_all_device<F,n_deriv>(
      weight_alg, func, sycl_data, tasks_device.begin(), tasks_device.end()
    );

    task_it = it;
  }

  // Receive XC terms from host
  if( not sycl_data.vxcinc_host )
    util::sycl_copy( nbf * nbf, VXC, sycl_data.vxc_device, syclQue, "VXC D2H" );

  util::sycl_copy( 1, EXC, sycl_data.exc_device, syclQue, "EXC D2H" );
  util::sycl_copy( 1, NEL, sycl_data.nel_device, syclQue, "NEL D2H" );

  // Symmetrize VXC
  for( int32_t j = 0;   j < nbf; ++j )
  for( int32_t i = j+1; i < nbf; ++i )
    VXC[ j + i*nbf ] = VXC[ i + j*nbf ];

}


#define SYCL_IMPL( F, ND ) \
template \
void process_batches_sycl_replicated_p<F, ND>(\
  XCWeightAlg            weight_alg,\
  const functional_type& func,\
  const BasisSet<F>&     basis,\
  const Molecule   &     mol,\
  const MolMeta    &     meta,\
  XCSyclData<F>    &     sycl_data,\
  std::vector< XCTask >& local_work,\
  const F*               P,\
  F*                     VXC,\
  F*                     exc,\
  F*                     n_el\
)

SYCL_IMPL( double, 0 );
SYCL_IMPL( double, 1 );

}
}
