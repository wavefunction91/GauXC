/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */

#include "incore_replicated_xc_device_integrator.hpp"
#include "device/local_device_work_driver.hpp"
#include "host/reference_local_host_work_driver.hpp"
#include <stdexcept>
#include "device/xc_device_aos_data.hpp"
#include <fstream>
#include <gauxc/util/unused.hpp>

#include "integrator_util/exx_screening.hpp"
#include "integrator_util/integral_bounds.hpp"

namespace GauXC  {
namespace detail {

template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
  eval_exx_grad_( int64_t m, int64_t n, const value_type* P,
             int64_t ldp, value_type* EXX_GRAD,
             const IntegratorSettingsEXX& settings ) {


  const auto& basis = this->load_balancer_->basis();

  // Check that P / K are sane
  const int64_t nbf = basis.nbf();
  if( m != n )
    GAUXC_GENERIC_EXCEPTION("P/K Must Be Square");
  if( m != nbf )
    GAUXC_GENERIC_EXCEPTION("P/K Must Have Same Dimension as Basis");
  if( ldp < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDP");

  // Allocate Device memory
  auto* lwd = dynamic_cast<LocalDeviceWorkDriver*>(this->local_work_driver_.get() );
  auto rt  = detail::as_device_runtime(this->load_balancer_->runtime());
  auto device_data_ptr = lwd->create_device_data(rt);

  std::vector<ValueType> exx_bfgrad(3*nbf, 0);

  GAUXC_MPI_CODE(MPI_Barrier(rt.comm());)

  this->timer_.time_op("XCIntegrator.EXX_GRAD_Screening", [&]() {
    exx_ek_screening_local_work_( basis, P, ldp, *device_data_ptr, settings);
  });


  // Get Tasks
  auto& tasks = this->load_balancer_->get_tasks();
  if( this->reduction_driver_->takes_device_memory() ) {

    // Compute local contributions to K and keep on device
    this->timer_.time_op("XCIntegrator.LocalWork_EXX_GRAD", [&](){
      exx_grad_local_work_( basis, P, ldp,
        tasks.begin(), tasks.end(), *device_data_ptr, settings);
      rt.device_backend()->master_queue_synchronize();
    });

    GAUXC_MPI_CODE(
    this->timer_.time_op("XCIntegrator.ImbalanceWait_EXX_GRAD",[&](){
      MPI_Barrier(rt.comm());
    });
    )

    // Reduce results in device memory
    this->timer_.time_op("XCIntegrator.Allreduce_EXX_GRAD", [&](){
      this->reduction_driver_->allreduce_inplace(
        device_data_ptr->exx_grad_device_data(), nbf, ReductionOp::Sum,
        device_data_ptr->queue());
    });

    // Receive K from host
    this->timer_.time_op("XCIntegrator.DeviceToHostCopy_EXX_GRAD",[&](){
      device_data_ptr->retrieve_exx_grad(exx_bfgrad.data());
    });

  } else {

    // Compute local contributions to K and retrieve
    // data from device
    this->timer_.time_op("XCIntegrator.LocalWork_EXX_GRAD", [&](){
      exx_grad_local_work_( basis, P, ldp,
        tasks.begin(), tasks.end(), *device_data_ptr, settings);
    });

    GAUXC_MPI_CODE(
    this->timer_.time_op("XCIntegrator.ImbalanceWait_EXX_GRAD",[&](){
      MPI_Barrier(rt.comm());
    });
    )

    this->timer_.time_op("XCIntegrator.DeviceToHostCopy_EXX_GRAD",[&](){
      device_data_ptr->retrieve_exx_grad(exx_bfgrad.data());
    });

    // Reduce Results in host mem
    this->timer_.time_op("XCIntegrator.Allreduce_EXX_GRAD", [&](){
      this->reduction_driver_->allreduce_inplace(exx_bfgrad.data(), 3*nbf, ReductionOp::Sum );
    });

  }

  // Sum gradient contribution of basis functions
  // to nuclei
  rt.device_backend()->master_queue_synchronize();
  auto& basis_map   = this->load_balancer_->basis_map();
for (size_t i = 0; i< basis.nshells(); i++) {
    const auto [b0, b1] = basis_map.shell_to_ao_range()[i];
    const auto iCenter = basis_map.shell_to_center()[i];
    for (size_t bf = b0; bf < b1; bf++) {
        // factor 2 is for bra-ket integral symmetry
        EXX_GRAD[3*iCenter  ] += 2*exx_bfgrad[bf  ];
        EXX_GRAD[3*iCenter+1] += 2*exx_bfgrad[bf+1*nbf];
        EXX_GRAD[3*iCenter+2] += 2*exx_bfgrad[bf+2*nbf];
    }

}
}

template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
  exx_grad_local_work_( const basis_type& basis, const value_type* P, int64_t ldp,
                       value_type* EXX_GRAD, int64_t nbf,
                       host_task_iterator task_begin, host_task_iterator task_end,
                       XCDeviceData& device_data,
                       const IntegratorSettingsEXX& settings ) {


  exx_local_work_(basis, P, ldp, task_begin, task_end, device_data, settings);
  auto rt  = detail::as_device_runtime(this->load_balancer_->runtime());
  rt.device_backend()->master_queue_synchronize();

  // Receive K from host
  this->timer_.time_op("XCIntegrator.DeviceToHostCopy_EXX_GRAD",[&](){
    device_data.retrieve_exx_integrands( EXX_GRAD, nbf );
  });

}

template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
  exx_grad_local_work_( const basis_type& basis, const value_type* P, int64_t ldp,
                       host_task_iterator task_begin, host_task_iterator task_end,
                       XCDeviceData& device_data,
                       const IntegratorSettingsEXX& settings ) {

  auto* lwd = dynamic_cast<LocalDeviceWorkDriver*>(this->local_work_driver_.get() );
  IntegratorSettingsSNLinK sn_link_settings;
  if( auto* tmp = dynamic_cast<const IntegratorSettingsSNLinK*>(&settings) ) {
    sn_link_settings = *tmp;
  }

  // Setup Aliases
  const auto nbf     = basis.nbf();
  const auto nshells = basis.nshells();


  // Get basis map and shell pairs
  auto& basis_map   = this->load_balancer_->basis_map();
  auto& shell_pairs = this->load_balancer_->shell_pairs();



  // Sort tasks
  auto task_comparator = []( const XCTask& a, const XCTask& b ) {
    return (a.points.size() * a.bfn_screening.nbe) > (b.points.size() * b.bfn_screening.nbe);
  };
  std::sort( task_begin, task_end, task_comparator );



  // Check that Partition Weights have been calculated
  auto& lb_state = this->load_balancer_->state();
  if( not lb_state.modified_weights_are_stored ) {
    GAUXC_GENERIC_EXCEPTION("Weights Have Not Been Modified");
  }

    task_end = std::stable_partition( task_begin, task_end,
      []( const auto& t ) { return t.cou_screening.shell_list.size() > 0; } );

  std::sort(task_begin,task_end,
    [](auto& a, auto& b){ return a.cou_screening.shell_pair_list.size() >
      b.cou_screening.shell_pair_list.size(); });

  // Populate submat maps
  device_data.populate_submat_maps( basis.nbf(), task_begin, task_end, basis_map );



  // Do EXX integration in task batches
  device_data.reset_allocations();
  device_data.allocate_static_data_exx_grad( nbf, nshells, shell_pairs.npairs(), shell_pairs.nprim_pair_total(), basis_map.max_l() );
  device_data.send_static_data_density_basis( P, ldp, nullptr, 0, nullptr, 0, nullptr, 0, basis );
  device_data.send_static_data_shell_pairs( basis, shell_pairs );

  // Zero integrands
  device_data.zero_exx_integrands();
  device_data.zero_exx_grad_integrands();

  // Processes batches in groups that saturadate available device memory
  integrator_term_tracker enabled_terms;
  enabled_terms.exx = true;
  enabled_terms.exx_grad = true;

  auto task_it = task_begin;
  while( task_it != task_end ) {

    // Determine next task batch, send relevant data to device (EXX only)
    task_it =
      device_data.generate_buffers( enabled_terms, basis_map, task_it, task_end );

#if 1
    /*** Process the batches ***/

    // Evaluate collocation gradient
    lwd->eval_collocation_gradient( &device_data );

    // Evaluate F(mu,i) = P(mu,nu) * B(nu,i)
    // mu runs over significant ek shells
    // nu runs over the bfn shell list
    // i runs over all points
    lwd->eval_exx_fmat( &device_data );

    // Compute G(mu,i) = w(i) * A(mu,nu,i) * F(nu,i)
    // mu/nu run over significant ek shells
    // i runs over all points
    lwd->eval_exx_gmat( &device_data, basis_map );

    // Increment dK(mu,nu)/dx += dB/dx(mu,i) * G(nu,i)
    // mu runs over bfn shell list
    // nu runs over ek shells
    // i runs over all points
    lwd->eval_exx_kgrad( &device_data );
#endif

  } // Loop over batches of batches

  // Contract derivative K matrices with density
  // to produce gradient:
  // bfgrad(mu)_x = sum_nu dK/dx(mu,nu) * DM(mu,nu)
  lwd->inc_exx_kgrad( &device_data );

}

}
}
