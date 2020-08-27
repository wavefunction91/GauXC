#include <gauxc/xc_integrator/xc_cuda_util.hpp>
#include <gauxc/util/cuda_util.hpp>

#include "cuda_weights.hpp"
#include "collocation_device.hpp"
//#include "cuda_zmat.hpp"
#include "integrator_common.hpp"
//#include "blas.hpp"
//#include "util.hpp"


namespace GauXC  {
namespace integrator::cuda {

using host_task_iterator = std::vector<XCTask>::iterator;

template <typename F>
using cuda_task_iterator = typename std::vector<XCTaskDevice<F>>::iterator;

template <typename F, size_t n_deriv>
void process_batches_cuda_replicated_all_device(
  XCWeightAlg            weight_alg,
  const functional_type& func,
  XCCudaData<F>&         cuda_data,
  cuda_task_iterator<F>  task_begin,
  cuda_task_iterator<F>  task_end
) {


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

  const size_t total_npts = 
    std::accumulate( task_begin, task_end, 0ul, 
                     []( const auto& a, const auto& b ) { return a + b.npts; } );
                      


  // Aliases
  cudaStream_t   master_stream = *cuda_data.master_stream;
  cublasHandle_t master_handle = *cuda_data.master_handle;
  magma_queue_t  master_queue  = *cuda_data.master_magma_queue;

  const auto*  rab_device          = cuda_data.rab_device;
  const auto*  coords_device       = cuda_data.coords_device;
  const auto*  points_device       = cuda_data.points_device_buffer;
  const auto*  iparent_device      = cuda_data.iparent_device_buffer;
  const auto*  dist_nearest_device = cuda_data.dist_nearest_buffer;

  auto* weights_device      = cuda_data.weights_device_buffer;
  auto* dist_scratch_device = cuda_data.dist_scratch_device;



  // Evaluate partition weights
  partition_weights_cuda_SoA( weight_alg, total_npts, cuda_data.natoms, 
                              points_device, iparent_device, dist_nearest_device,
                              rab_device, coords_device, weights_device, 
                              dist_scratch_device, master_stream );


}


template <typename F, size_t n_deriv>
void process_batches_cuda_replicated_p(
  XCWeightAlg            weight_alg,
  const functional_type& func,
  const BasisSet<F>&     basis,
  const Molecule   &     mol,
  const MolMeta    &     meta,
  XCCudaData<F>    &     cuda_data,
  std::vector< XCTask >& tasks,
  const F*               P,
  F*                     VXC,
  F*                     EXC,
  F*                     NEL
) {


  auto task_comparator = []( const XCTask& a, const XCTask& b ) {
    return (a.points.size() * a.nbe) > (b.points.size() * b.nbe);
  };
  std::sort( tasks.begin(), tasks.end(), task_comparator );


  const auto nbf    = basis.nbf();
  const auto natoms = meta.natoms();
  // Send static data to the device

  // Density
  if( not cuda_data.denpack_host )
    util::cuda_copy( nbf * nbf, cuda_data.dmat_device, P, "P H2D" );

  // Shells: TODO avoid copy?
  std::vector<Shell<F>> shells( basis );
  util::cuda_copy( shells.size(), cuda_data.shells_device, shells.data(),
                   "Shells H2D" );

  // RAB
  util::cuda_copy( natoms * natoms, cuda_data.rab_device, meta.rab().data(),
                   "RAB H2D" );

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
  auto task_it = tasks.begin();
  while( task_it != tasks.end() ) {

    // Determine next task batch, send relevant data to device
    auto [it, tasks_device] = 
      cuda_data.generate_buffers( basis, task_it, tasks.end() );

    // TODO Process the batches
    process_batches_cuda_replicated_all_device<F,n_deriv>( 
      weight_alg, func, cuda_data, tasks_device.begin(), tasks_device.end() 
    );

    task_it = it;

  }

  // Receive XC terms from host
  if( not cuda_data.vxcinc_host )
    util::cuda_copy( nbf * nbf, VXC, cuda_data.vxc_device, "VXC D2H" );

  util::cuda_copy( 1, EXC, cuda_data.exc_device, "EXC D2H" );
  util::cuda_copy( 1, NEL, cuda_data.nel_device, "NEL D2H" );

}


#define CUDA_IMPL( F, ND ) \
template \
void process_batches_cuda_replicated_p<F, ND>(\
  XCWeightAlg            weight_alg,\
  const functional_type& func,\
  const BasisSet<F>&     basis,\
  const Molecule   &     mol,\
  const MolMeta    &     meta,\
  XCCudaData<F>    &     cuda_data,\
  std::vector< XCTask >& local_work,\
  const F*               P,\
  F*                     VXC,\
  F*                     exc,\
  F*                     n_el\
) 

CUDA_IMPL( double, 0 );
CUDA_IMPL( double, 1 );

}
}

