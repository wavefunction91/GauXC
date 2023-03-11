/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for retails
 */
#pragma once

#include <gauxc/load_balancer.hpp>

namespace GauXC  {
namespace detail {

class LoadBalancerImpl {

public:

  using basis_type = BasisSet<double>;

protected:

  RuntimeEnvironment          runtime_;
  std::shared_ptr<Molecule>   mol_;
  std::shared_ptr<MolGrid>    mg_;
  std::shared_ptr<basis_type> basis_;
  std::shared_ptr<MolMeta>    molmeta_;

  std::vector< XCTask >     local_tasks_;

  LoadBalancerState         state_;

  util::Timer               timer_;

  size_t                    pad_value_;

  virtual std::vector< XCTask > create_local_tasks_() const = 0;

public:

  LoadBalancerImpl() = delete;

  LoadBalancerImpl( const RuntimeEnvironment&, const Molecule&, const MolGrid& mg,  
    const basis_type&, size_t pv );
  LoadBalancerImpl( const RuntimeEnvironment&, const Molecule&, const MolGrid& mg,  
    const basis_type&, const MolMeta&, size_t pv );
  LoadBalancerImpl( const RuntimeEnvironment&, const Molecule&, const MolGrid& mg,  
    const basis_type&, std::shared_ptr<MolMeta>, size_t pv );

  LoadBalancerImpl( const LoadBalancerImpl& );
  LoadBalancerImpl( LoadBalancerImpl&& ) noexcept;

  virtual ~LoadBalancerImpl() noexcept;

  const std::vector< XCTask >& get_tasks() const;
        std::vector< XCTask >& get_tasks()      ;

  void rebalance_weights();
  void rebalance_exc_vxc();
  void rebalance_exx();

  const util::Timer& get_timings() const;

  size_t max_npts()       const;
  size_t max_nbe()        const;
  size_t max_npts_x_nbe() const;
  size_t pad_value()      const;

  const Molecule& molecule() const;
  const MolMeta&  molmeta()  const;
  const basis_type& basis()  const;
  const RuntimeEnvironment& runtime() const;

  LoadBalancerState& state();

  virtual std::unique_ptr<LoadBalancerImpl> clone() const = 0;

};

}
}
