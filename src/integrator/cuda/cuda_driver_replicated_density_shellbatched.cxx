#include <gauxc/xc_integrator/xc_cuda_util.hpp>
#include <gauxc/util/cuda_util.hpp>
#include <gauxc/util/unused.hpp>

#include "cuda/cuda_weights.hpp"
#include "cuda/collocation_device.hpp"
#include "cuda/cuda_pack_density.hpp"
#include "cuda/cuda_inc_potential.hpp"
#include "cuda/cuda_eval_denvars.hpp"
#include "cuda/cuda_zmat.hpp"
#include "integrator_common.hpp"
  
#include "cuda/cublas_extensions.hpp"

#include "host/util.hpp"

namespace GauXC  {
namespace integrator::cuda {

using namespace GauXC::cuda::blas;


template <typename F, size_t n_deriv>
void process_batches_cuda_replicated_density_shellbatched_p(
  XCWeightAlg            weight_alg,
  const functional_type& func,
  const BasisSet<F>&     basis,
  const Molecule   &     mol,
  const MolMeta    &     meta,
  XCCudaData<F>    &     cuda_data,
  host_task_iterator     local_work_begin,
  host_task_iterator     local_work_end,
  const F*               P,
  F*                     VXC,
  F*                     EXC,
  F*                     NEL
) {

  std::cout << "IN SHELL BATCHED\n" << std::flush;


  // Zero out final results
  *EXC = 0.;
  *NEL = 0.;
  for( auto i = 0; i < basis.nbf() * basis.nbf(); ++i )
    VXC[i] = 0.;

#if 0
  size_t nbf     = basis.nbf();
  size_t nshells = basis.nshells();
  size_t natoms  = mol.size();

  // Allocate static quantities on device stack
  cuda_data.allocate_static_data( natoms, n_deriv, nbf, nshells );

  process_batches_cuda_replicated_density_incore_p<F,n_deriv>(
    weight_alg, func, basis, mol, meta, cuda_data, 
    local_work_begin, local_work_end, P, VXC, EXC, NEL
  );
#else

  auto nbe_comparator = []( const auto& task_a, const auto& task_b ) {
    return task_a.nbe < task_b.nbe;
  };

  auto task_begin = local_work_begin;
  while( task_begin != local_work_end ) {

    // Find task with largest NBE
    auto max_task = std::max_element( task_begin, local_work_end, nbe_comparator );
    const auto   max_shell_list = max_task->shell_list; // copy for reset
    const size_t nbe     = max_task->nbe;
    const size_t nshells = max_shell_list.size();
    const size_t natoms  = mol.size();

    // Partition tasks into those which are subsets of max_task and those 
    // that aren't
    auto subset_of_max = [&]( const auto& t ) {
      return std::includes( max_shell_list.begin(),
                            max_shell_list.end(),
                            t.shell_list.begin(),
                            t.shell_list.end() );
    };

    auto task_end = std::partition( task_begin, local_work_end, subset_of_max );


    // Create subobjects for batch processing

    std::cout << "MAX TASK HAS:"   << std::endl
              << "  NSHELLS    = " << nshells << std::endl
              << "  NBE        = " << nbe     << std::endl;
    std::cout << "FOUND " << std::distance( task_begin, task_end ) << " SUBTASKS " << std::endl;
              

    // Extract subbasis
    BasisSet<F> basis_subset; basis_subset.reserve(nshells);
    for( auto i : max_shell_list ) {
      basis_subset.emplace_back( basis.at(i) );
    }
    basis_subset.generate_shell_to_ao();

    // XXX: DEBUG
    //std::cout << "SHELL LISTS BEFORE RECALC" << std::endl;
    //for( auto _it = task_begin; _it != task_end; ++_it ) {
    //  std::cout << "TASK " << std::distance( task_begin, _it ) << ": ";
    //  for( auto j : _it->shell_list ) std::cout << j << ", ";
    //  std::cout << std::endl;
    //}


    // Recalculate shell_list based on subbasis
    for( auto i = 0ul; i < nshells; ++i ) {
      const auto shell_idx = max_shell_list[i];
      for( auto _it = task_begin; _it != task_end; ++_it )
      for( auto j = 0ul; j < _it->shell_list.size(); ++j ) 
      if( _it->shell_list[j] == shell_idx ) {
        _it->shell_list[j] = i;
      }
    }
    
    // XXX: DEBUG
    //std::cout << "SHELL LISTS AFTER RECALC" << std::endl;
    //for( auto _it = task_begin; _it != task_end; ++_it ) {
    //  std::cout << "TASK " << std::distance( task_begin, _it ) << ": ";
    //  for( auto j : _it->shell_list ) std::cout << j << ", ";
    //  std::cout << std::endl;
    //}

#if 1

    // Allocate temporary submatrices
    std::vector<F> P_submat_host(nbe*nbe), VXC_submat_host(nbe*nbe);
    F EXC_tmp, NEL_tmp;
    F* P_submat   = P_submat_host.data();
    F* VXC_submat = VXC_submat_host.data();

    // Extract subdensity
    auto [max_submat_cut, foo] = 
      integrator::gen_compressed_submat_map( basis, max_shell_list, 
        basis.nbf(), basis.nbf() );

    detail::submat_set( basis.nbf(), basis.nbf(), nbe, nbe, P, basis.nbf(), 
                        P_submat, nbe, max_submat_cut );
   

    // Allocate static quantities on device stack
    cuda_data.allocate_static_data( natoms, n_deriv, nbe, nshells );

    // Process batches on device with subobjects
    process_batches_cuda_replicated_density_incore_p<F,n_deriv>(
      weight_alg, func, basis_subset, mol, meta, cuda_data, 
      task_begin, task_end, P_submat, VXC_submat, &EXC_tmp, &NEL_tmp
    );

    // Update full quantities
    *EXC += EXC_tmp;
    *NEL += NEL_tmp;
    detail::inc_by_submat( basis.nbf(), basis.nbf(), nbe, nbe, VXC, basis.nbf(), 
                           VXC_submat, nbe, max_submat_cut );

#endif

    // Reset shell_list to be wrt full basis
    for( auto i = 0ul; i < nshells; ++i ) {
      const auto shell_idx = max_shell_list[i];
      for( auto _it = task_begin; _it != task_end; ++_it )
      for( auto j = 0ul; j < _it->shell_list.size(); ++j ) 
      if( _it->shell_list[j] == i ) {
        _it->shell_list[j] = shell_idx;
      }
    }

    // XXX: DEBUG
    //std::cout << "SHELL LISTS AFTER RESET" << std::endl;
    //for( auto _it = task_begin; _it != task_end; ++_it ) {
    //  std::cout << "TASK " << std::distance( task_begin, _it ) << ": ";
    //  for( auto j : _it->shell_list ) std::cout << j << ", ";
    //  std::cout << std::endl;
    //}

    // Update task iterator for next set of batches
    task_begin = task_end;

  }


#endif

}


#define CUDA_IMPL( F, ND ) \
template \
void process_batches_cuda_replicated_density_shellbatched_p<F, ND>(\
  XCWeightAlg            weight_alg,\
  const functional_type& func,\
  const BasisSet<F>&     basis,\
  const Molecule   &     mol,\
  const MolMeta    &     meta,\
  XCCudaData<F>    &     cuda_data,\
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

