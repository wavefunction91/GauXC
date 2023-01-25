/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for retails
 */
#pragma once
#include <gauxc/load_balancer.hpp>
#include <gauxc/util/timer.hpp>

namespace GauXC {

namespace detail {
  class MolecularWeightsImpl;

}

struct MolecularWeightsSettings { };

class MolecularWeights {

public:

  using load_balancer_type = LoadBalancer;
  using load_balancer_reference = load_balancer_type&;

private:

  using pimpl_type = detail::MolecularWeightsImpl;
  using pimpl_ptr_type = std::unique_ptr<pimpl_type>;
  pimpl_ptr_type pimpl_;

public:

  MolecularWeights() = delete;
  ~MolecularWeights() noexcept;

  MolecularWeights( pimpl_ptr_type&& ptr );

  MolecularWeights( const MolecularWeights& ) = delete;
  MolecularWeights( MolecularWeights&& ) noexcept;

  void modify_weights(load_balancer_reference lb) const;
  const util::Timer& get_timings() const;

};



class MolecularWeightsFactory {

public:

    MolecularWeightsFactory() = delete;
    MolecularWeightsFactory( ExecutionSpace ex, 
                             std::string local_work_kernel_name,
                             MolecularWeightsSettings s);


    std::shared_ptr<MolecularWeights> get_shared_instance();

    inline MolecularWeights get_instance(){
      return MolecularWeights( std::move( *get_shared_instance() ) );
    };

private:

    ExecutionSpace ex_;
    std::string lwd_kernel_;
    MolecularWeightsSettings settings_;

};

}
