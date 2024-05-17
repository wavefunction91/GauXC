/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once

#include <gauxc/xc_integrator/xc_integrator_impl.hpp>

namespace GauXC  {
namespace detail {

template <typename ValueType>
class ReplicatedXCIntegratorImpl;


/** XCIntegrator implementation for replicated inputs
 *
 *  Expects for the passed MatrixType to be convertable to a 
 *  dense matrix and that the inputs are replicted on every rank.
 */
template <typename MatrixType>
class ReplicatedXCIntegrator : public XCIntegratorImpl<MatrixType> {

public:

  using matrix_type    = typename XCIntegratorImpl<MatrixType>::matrix_type;
  using value_type     = typename XCIntegratorImpl<MatrixType>::value_type;
  using exc_vxc_type_rks   = typename XCIntegratorImpl<MatrixType>::exc_vxc_type_rks;
  using exc_vxc_type_uks   = typename XCIntegratorImpl<MatrixType>::exc_vxc_type_uks;
  using exc_vxc_type_gks   = typename XCIntegratorImpl<MatrixType>::exc_vxc_type_gks;
  using exc_grad_type  = typename XCIntegratorImpl<MatrixType>::exc_grad_type;
  using exx_type       = typename XCIntegratorImpl<MatrixType>::exx_type;
  using atomic_overlap_type   = typename XCIntegratorImpl<MatrixType>::atomic_overlap_type;

private:

  using pimpl_type = ReplicatedXCIntegratorImpl<value_type>;
  std::unique_ptr< pimpl_type > pimpl_;

  value_type    integrate_den_( const MatrixType& ) override;
  exc_vxc_type_rks  eval_exc_vxc_ ( const MatrixType&, const IntegratorSettingsXC& ) override;
  exc_vxc_type_uks  eval_exc_vxc_ ( const MatrixType&, const MatrixType&, const IntegratorSettingsXC&) override;
  exc_vxc_type_gks  eval_exc_vxc_ ( const MatrixType&, const MatrixType&, const MatrixType&, const MatrixType&, const IntegratorSettingsXC& ) override;
  exc_grad_type eval_exc_grad_( const MatrixType& ) override;
  exx_type      eval_exx_     ( const MatrixType&, const IntegratorSettingsEXX& ) override;
  atomic_overlap_type eval_atomic_overlap_(int64_t iAtom) override;

  const util::Timer& get_timings_() const override;
  const LoadBalancer& get_load_balancer_() const override;
  LoadBalancer& get_load_balancer_() override;

public:

  ReplicatedXCIntegrator();
  ReplicatedXCIntegrator( std::unique_ptr<pimpl_type>&& );

  ~ReplicatedXCIntegrator() noexcept;

  ReplicatedXCIntegrator( const ReplicatedXCIntegrator& ) = delete;
  ReplicatedXCIntegrator( ReplicatedXCIntegrator&& ) noexcept;

};


}
}
