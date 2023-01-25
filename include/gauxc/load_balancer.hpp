/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for retails
 */
#pragma once

#include <gauxc/molgrid.hpp>
#include <gauxc/molmeta.hpp>
#include <gauxc/basisset.hpp>
#include <gauxc/xc_task.hpp>
#include <gauxc/util/timer.hpp>
#include <gauxc/runtime_environment.hpp>

namespace GauXC {

namespace detail {
  class LoadBalancerImpl;
}

struct LoadBalancerState {
  bool modified_weights_are_stored = false;
};

class LoadBalancer {

  using basis_type = BasisSet<double>;
  using pimpl_type = detail::LoadBalancerImpl;
  std::unique_ptr<pimpl_type> pimpl_;

public:

  LoadBalancer();

  LoadBalancer( std::unique_ptr<pimpl_type>&& pimpl );

  LoadBalancer( const LoadBalancer& );
  LoadBalancer( LoadBalancer&& ) noexcept;

  ~LoadBalancer() noexcept;

  const std::vector<XCTask>& get_tasks() const;
        std::vector<XCTask>& get_tasks()      ;
  
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

  bool operator==( const LoadBalancer& ) const;
};



class LoadBalancerFactory {

public:

  LoadBalancerFactory() = delete;

  LoadBalancerFactory( ExecutionSpace ex, std::string kernel_name );

  LoadBalancer get_instance( const RuntimeEnvironment& rt, 
    const Molecule& mol, const MolGrid& mg, const BasisSet<double>&,
    size_t pad_val = 1 );
  std::shared_ptr<LoadBalancer> get_shared_instance( 
    const RuntimeEnvironment& rt,
    const Molecule& mol, const MolGrid& mg, const BasisSet<double>&,
    size_t pad_val = 1 );

private:

  ExecutionSpace ex_;
  std::string    kernel_name_;

};


}
