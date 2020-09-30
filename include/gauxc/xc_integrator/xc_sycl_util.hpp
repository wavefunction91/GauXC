#pragma once
#include <gauxc/xc_integrator/xc_sycl_data.hpp>

#include <gauxc/xc_integrator.hpp>

#ifdef GAUXC_ENABLE_SYCL
namespace GauXC  {
namespace integrator {
namespace sycl {


template <typename F, size_t n_deriv>
void process_batches_sycl_replicated_p(
  XCWeightAlg            weight_alg,
  const functional_type& func,
  const BasisSet<F>&     basis,
  const Molecule   &     mol,
  const MolMeta    &     meta,
  XCSyclData<F>    &     sycl_data,
  std::vector< XCTask >& local_work,
  const F*               P,
  F*                     VXC,
  F*                     exc,
  F*                     n_el
);


template <typename F, typename... Args>
inline void process_batches_sycl_replicated_p( size_t n_deriv, Args&&... args ) {
  if( n_deriv == 0 )
    process_batches_sycl_replicated_p<F,0>( std::forward<Args>(args)... );
  else if( n_deriv == 1 )
    process_batches_sycl_replicated_p<F,1>( std::forward<Args>(args)... );
  else
    throw std::runtime_error("MGGA NYI");
}

}
}
}
#endif
