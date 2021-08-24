#include <memory>
#include <gauxc/util/cuda_util.hpp>
#include <gauxc/util/unused.hpp>

#include "device/cuda/cuda_weights.hpp"
#include "device/cuda/collocation_device.hpp"
#include "device/cuda/cuda_pack_density.hpp"
#include "device/cuda/cuda_inc_potential.hpp"
#include "device/cuda/cuda_eval_denvars.hpp"
#include "device/cuda/cuda_zmat.hpp"
#include "common/integrator_common.hpp"
  
#include "device/cuda/cublas_extensions.hpp"
#include "device/cuda/local_work_replicated_incore_exc_vxc.hpp"

#include "device/cuda/xc_cuda_data.hpp"

namespace GauXC  {

namespace integrator::cuda {

using namespace GauXC::cuda::blas;


template <typename F>
using cuda_task_iterator = typename std::vector<XCTaskDevice<F>>::iterator;

template <typename F, size_t n_deriv>
void local_work_replicated_density_incore_exc_vxc(
  XCWeightAlg            weight_alg,
  const functional_type& func,
  XCCudaData<F>&         cuda_data,
  cuda_task_iterator<F>  task_begin,
  cuda_task_iterator<F>  task_end
) {

  const auto ntasks = std::distance( task_begin, task_end );
  const auto nbf    = cuda_data.nbf;

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

  util::unused( min_nbe, min_npts, min_nshells );

  const size_t total_npts = 
    std::accumulate( task_begin, task_end, 0ul, 
                     []( const auto& a, const auto& b ) { return a + b.npts; } );


  // Aliases
  cudaStream_t   master_stream = *cuda_data.master_stream;
  cublasHandle_t master_handle = *cuda_data.master_handle;

#ifdef GAUXC_ENABLE_MAGMA
  magma_queue_t  master_queue  = *cuda_data.master_magma_queue;
#endif

  auto* dmat_device         = cuda_data.dmat_device;

  auto* shells_device       = cuda_data.shells_device;
  auto* tasks_device        = cuda_data.device_tasks;
  auto* dmat_array_device   = cuda_data.dmat_array_device;
  auto* zmat_array_device   = cuda_data.zmat_array_device;
  auto* bf_array_device     = cuda_data.bf_array_device;
  auto* weights_device      = cuda_data.weights_device_buffer;
  auto* dist_scratch_device = cuda_data.dist_scratch_device;

  auto* den_eval_device     = cuda_data.den_eval_device;
  auto* dden_x_eval_device  = cuda_data.den_x_eval_device;
  auto* dden_y_eval_device  = cuda_data.den_y_eval_device;
  auto* dden_z_eval_device  = cuda_data.den_z_eval_device;

  auto* eps_eval_device     = cuda_data.eps_eval_device;
  auto* gamma_eval_device   = cuda_data.gamma_eval_device;
  auto* vrho_eval_device    = cuda_data.vrho_eval_device;
  auto* vgamma_eval_device  = cuda_data.vgamma_eval_device;


  auto* exc_device     = cuda_data.exc_device;
  auto* vxc_device     = cuda_data.vxc_device;
  auto* nel_device     = cuda_data.nel_device;
  auto* acc_scr_device = cuda_data.acc_scr_device;

  auto* m_array_device      = cuda_data.m_array_device;
  auto* n_array_device      = cuda_data.n_array_device;
  auto* k_array_device      = cuda_data.k_array_device;
  auto* lda_array_device    = cuda_data.lda_array_device;
  auto* ldb_array_device    = cuda_data.ldb_array_device;
  auto* ldc_array_device    = cuda_data.ldc_array_device;


  const auto* rab_device          = cuda_data.rab_device;
  const auto* coords_device       = cuda_data.coords_device;
  const auto* points_device       = cuda_data.points_device_buffer;
  const auto* iparent_device      = cuda_data.iparent_device_buffer;
  const auto* dist_nearest_device = cuda_data.dist_nearest_buffer;




  // Evaluate Partition Weights
  partition_weights_cuda_SoA( weight_alg, total_npts, cuda_data.LDatoms, cuda_data.natoms, 
                              points_device, iparent_device, dist_nearest_device,
                              rab_device, coords_device, weights_device, 
                              dist_scratch_device, master_stream );


  // Evaluate Collocation
  if constexpr ( n_deriv == 1 )
    eval_collocation_masked_combined_deriv1( ntasks, max_npts, max_nshells,
                                             shells_device, tasks_device,
                                             master_stream );
  else
    eval_collocation_masked_combined( ntasks, max_npts, max_nshells, shells_device, 
                                      tasks_device, master_stream );

  // Pack Density Submatrices
  task_pack_density_matrix( ntasks, tasks_device, dmat_device, nbf, master_stream );


  // Form Z = P * X
  if( cuda_data.batch_l3_blas ) {

#ifdef GAUXC_ENABLE_MAGMA

    magmablas_dgemm_vbatched( MagmaNoTrans, MagmaNoTrans,
                              m_array_device, n_array_device, k_array_device,
                              1., bf_array_device, ldb_array_device,
                              dmat_array_device, lda_array_device,
                              0., zmat_array_device, ldc_array_device,
                              ntasks, master_queue );

#else

    throw std::runtime_error("BATCHED BLAS API NOT SUPPORTED");

#endif

  } else {

    int nstream = cuda_data.blas_streams.size();

    // Wait for collocation etc
    util::cuda_event master_event;
    master_event.record( master_stream );
    for( int iS = 0; iS < nstream; ++iS ) 
      cuda_data.blas_streams[iS].wait( master_event );

    // Do GEMM in round-robin
    for( auto iT = 0; iT < ntasks; ++iT ) {
      auto& task = *(task_begin + iT);
      gemm( cuda_data.blas_handles[iT % nstream], CUBLAS_OP_N, CUBLAS_OP_N,
            task.npts, task.nbe, task.nbe, 1., task.bf, task.npts,
            task.nbe_scr, task.nbe, 0., task.zmat, task.npts );
    }

    // Record completion of BLAS ops
    std::vector< util::cuda_event > blas_events( nstream );
    for( int iS = 0; iS < nstream; ++iS ) 
      blas_events[iS].record( cuda_data.blas_streams[iS] );

    // Wait on master stream for all BLAS ops to complete
    for( int iS = 0; iS < nstream; ++iS )
      cuda_data.master_stream->wait( blas_events[iS] );

  }
                

  
  // Zero UVars
  util::cuda_set_zero_async( total_npts, den_eval_device, master_stream, "DenZero" );
  if( func.is_gga() ) {
    util::cuda_set_zero_async( total_npts, dden_x_eval_device, master_stream, 
                               "DenXZero" );
    util::cuda_set_zero_async( total_npts, dden_y_eval_device, master_stream, 
                               "DenYZero" );
    util::cuda_set_zero_async( total_npts, dden_z_eval_device, master_stream, 
                               "DenZZero" );
  }

  // Evaluate UVars
  if( func.is_gga() ) {
    eval_uvars_gga_device( ntasks, max_nbe, max_npts, tasks_device, master_stream );
    eval_vvars_gga_device( total_npts, dden_x_eval_device, dden_y_eval_device,
                           dden_z_eval_device, gamma_eval_device, master_stream );
  } else {
    eval_uvars_lda_device( ntasks, max_nbe, max_npts, tasks_device, master_stream );
  }

  // Evaluate XC Functional
  if( func.is_gga() )
    func.eval_exc_vxc_device( total_npts, den_eval_device, gamma_eval_device, 
                              eps_eval_device, vrho_eval_device, 
                              vgamma_eval_device, master_stream );
  else
    func.eval_exc_vxc_device( total_npts, den_eval_device, eps_eval_device, 
                              vrho_eval_device, master_stream );


  // Factor weights into XC output
  hadamard_product( master_handle, total_npts, 1, weights_device, 1,
                    eps_eval_device, 1 );
  hadamard_product( master_handle, total_npts, 1, weights_device, 1,
                    vrho_eval_device, 1 );
  if( func.is_gga() ) 
    hadamard_product( master_handle, total_npts, 1, weights_device, 1,
                      vgamma_eval_device, 1 );

  // Accumulate EXC / NEL
  gdot( master_handle, total_npts, weights_device, 1,
        den_eval_device, 1, acc_scr_device, nel_device );
  gdot( master_handle, total_npts, eps_eval_device, 1,
        den_eval_device, 1, acc_scr_device, exc_device );
      
  // Evaluate Z Matrix
  if( func.is_gga() )
    zmat_gga_cuda( ntasks, max_nbe, max_npts, tasks_device, master_stream );
  else
    zmat_lda_cuda( ntasks, max_nbe, max_npts, tasks_device, master_stream );
  


  // Accumulate packed VXC = X * Z**T + Z * X**T

  
  if( cuda_data.batch_l3_blas ) {

#ifdef GAUXC_ENABLE_MAGMA

    // XXX: Only updates LT
    magmablas_dsyr2k_vbatched( MagmaLower, MagmaTrans, 
                               n_array_device, m_array_device,
                               1., bf_array_device, ldb_array_device,
                               zmat_array_device, ldc_array_device,
                               0., dmat_array_device, lda_array_device,
                               ntasks, master_queue );

#else

    throw std::runtime_error("BATCHED BLAS API NOT SUPPORTED");

#endif
  } else {

    int nstream = cuda_data.blas_streams.size();

    // Wait for zmat, etc
    util::cuda_event master_event;
    master_event.record( master_stream );
    for( int iS = 0; iS < nstream; ++iS ) 
      cuda_data.blas_streams[iS].wait( master_event );

    // Do SYR2K in round-robin
    for( auto iT = 0; iT < ntasks; ++iT ) {
      auto& task = *(task_begin + iT);
      syr2k( cuda_data.blas_handles[iT % nstream], CUBLAS_FILL_MODE_LOWER, 
             CUBLAS_OP_T, task.nbe, task.npts, 1., task.bf, task.npts,
             task.zmat, task.npts, 0., task.nbe_scr, task.nbe );
    }

    // Record completion of BLAS ops
    std::vector< util::cuda_event > blas_events( nstream );
    for( int iS = 0; iS < nstream; ++iS ) 
      blas_events[iS].record( cuda_data.blas_streams[iS] );

    // Wait on master stream for all BLAS ops to complete
    for( int iS = 0; iS < nstream; ++iS )
      cuda_data.master_stream->wait( blas_events[iS] );
  }

  // Increment global VXC
  task_inc_potential( ntasks, tasks_device, vxc_device, nbf, master_stream );


  // Synchronize on master stream
  // XXX: There's no lifetime issues in this driver, should look into
  //      avoid this sync to allow for overlap with the host packing 
  cudaStreamSynchronize( master_stream );

}


template <typename F, size_t n_deriv>
void local_work_replicated_incore_exc_vxc_impl(
  XCWeightAlg            weight_alg,
  XCIntegratorState      state,
  const functional_type& func,
  const BasisSet<F>&     basis,
  const Molecule   &     mol,
  const MolMeta    &     meta,
  XCDeviceData<F>  &     device_data,
  host_task_iterator     local_work_begin,
  host_task_iterator     local_work_end,
  const F*               P,
  F*                     VXC,
  F*                     EXC,
  F*                     NEL
) {

  auto& cuda_data = dynamic_cast< XCCudaData<F>& >( device_data );

  auto task_comparator = []( const XCTask& a, const XCTask& b ) {
    return (a.points.size() * a.nbe) > (b.points.size() * b.nbe);
  };
  std::sort( local_work_begin, local_work_end, task_comparator );


  const auto nbf     = basis.nbf();
  const auto natoms  = meta.natoms();
  const auto LDatoms = cuda_data.LDatoms;

  // Send static data to the device

  // Density
  util::cuda_copy( nbf * nbf, cuda_data.dmat_device, P, "P H2D" );

  // Shells: TODO avoid host copy?
  std::vector<Shell<F>> shells( basis );
  util::cuda_copy( shells.size(), cuda_data.shells_device, shells.data(),
                   "Shells H2D" );

  // RAB
  util::cuda_copy_2d( cuda_data.rab_device, LDatoms * sizeof(F),
                      meta.rab().data(), natoms * sizeof(F),
                      natoms * sizeof(F), natoms, "RAB H2D");
  // This could probably happen on the host
  cuda_reciprocal(natoms * LDatoms, cuda_data.rab_device, 0);

  // Atomic coordinates 
  std::vector<double> coords( 3*natoms );
  for( auto i = 0ul; i < natoms; ++i ) {
    coords[ 3*i + 0 ] = mol[i].x;
    coords[ 3*i + 1 ] = mol[i].y;
    coords[ 3*i + 2 ] = mol[i].z;
  }
  util::cuda_copy( 3 * natoms, cuda_data.coords_device, coords.data(),
                   "Coords H2D" );


  // Zero out XC quantities
  util::cuda_set_zero( nbf * nbf, cuda_data.vxc_device, "VXC Zero" ); 
  util::cuda_set_zero( 1        , cuda_data.exc_device, "EXC Zero" ); 
  util::cuda_set_zero( 1        , cuda_data.nel_device, "NEL Zero" ); 



  // Processes batches in groups that saturadate available device memory
  auto task_it = local_work_begin;
  while( task_it != local_work_end ) {

    // Determine next task batch, send relevant data to device
    auto [it, tasks_device] = 
      cuda_data.generate_buffers( basis, task_it, local_work_end );


    // Process the batches
    local_work_replicated_density_incore_exc_vxc<F,n_deriv>( 
      weight_alg, func, cuda_data, tasks_device.begin(), tasks_device.end() 
    );

    task_it = it;

  }

  // Receive XC terms from host
  util::cuda_copy( nbf * nbf, VXC, cuda_data.vxc_device, "VXC D2H" );

  util::cuda_copy( 1, EXC, cuda_data.exc_device, "EXC D2H" );
  util::cuda_copy( 1, NEL, cuda_data.nel_device, "NEL D2H" );

  // Symmetrize VXC
  for( int32_t j = 0;   j < nbf; ++j )
  for( int32_t i = j+1; i < nbf; ++i )
    VXC[ j + i*nbf ] = VXC[ i + j*nbf ];

}


#define CUDA_IMPL( F, ND ) \
template \
void local_work_replicated_incore_exc_vxc_impl<F,ND>(\
  XCWeightAlg            weight_alg,\
  XCIntegratorState      state,\
  const functional_type& func,\
  const BasisSet<F>&     basis,\
  const Molecule   &     mol,\
  const MolMeta    &     meta,\
  XCDeviceData<F>  &     device_data,\
  host_task_iterator     local_work_begin,\
  host_task_iterator     local_work_end,\
  const F*               P,\
  F*                     VXC,\
  F*                     exc,\
  F*                     n_el\
) 

CUDA_IMPL( double, 0 );
CUDA_IMPL( double, 1 );

}
}


