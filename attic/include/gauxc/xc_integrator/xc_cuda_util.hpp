#pragma once
#include <gauxc/xc_integrator/xc_cuda_data.hpp>
#include <gauxc/xc_integrator.hpp>

#include <gauxc/util/timer.hpp>

#ifdef GAUXC_ENABLE_CUDA
namespace GauXC  {
namespace integrator {
namespace cuda {

using host_task_iterator = std::vector<XCTask>::iterator;

template <typename F, size_t n_deriv>
void process_batches_cuda_replicated_density_incore_p(
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
  F*                     exc,
  F*                     n_el
);


template <typename F, typename... Args>
inline void process_batches_cuda_replicated_density_incore_p( size_t n_deriv, Args&&... args ) {
  if( n_deriv == 0 )
    process_batches_cuda_replicated_density_incore_p<F,0>( std::forward<Args>(args)... );
  else if( n_deriv == 1 )
    process_batches_cuda_replicated_density_incore_p<F,1>( std::forward<Args>(args)... );
  else
    throw std::runtime_error("MGGA NYI");
}



template <typename F, size_t n_deriv>
void process_batches_cuda_replicated_density_shellbatched_p(
  util::Timer&           timer,
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
  F*                     exc,
  F*                     n_el
);

#ifdef GAUXC_ENABLE_NCCL
template <typename F>
void device_allreduce(
  ncclComm_t             nccl_comm,
  XCCudaData<F>    &     cuda_data
); 
#endif

template <typename F>
void device_transfer(
  XCCudaData<F>    &     cuda_data,
  F*                     VXC,
  F*                     EXC,
  F*                     NEL
); 

template <typename F, typename... Args>
inline void process_batches_cuda_replicated_density_shellbatched_p( size_t n_deriv, Args&&... args ) {
  if( n_deriv == 0 )
    process_batches_cuda_replicated_density_shellbatched_p<F,0>( std::forward<Args>(args)... );
  else if( n_deriv == 1 )
    process_batches_cuda_replicated_density_shellbatched_p<F,1>( std::forward<Args>(args)... );
  else
    throw std::runtime_error("MGGA NYI");
}

}
}
}
#endif
