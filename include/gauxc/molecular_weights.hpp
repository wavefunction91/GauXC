/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy).
 *
 * (c) 2024-2025, Microsoft Corporation
 *
 * All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include <gauxc/load_balancer.hpp>
#include <gauxc/util/timer.hpp>
#include <gauxc/enums.hpp>

namespace GauXC {

namespace detail {
  /// Forward declaration of MolecularWeights implementation.
  class MolecularWeightsImpl;
}

/**
 *  @brief Settings for molecular weight partitioning.
 *
 *  Controls the weight partitioning algorithm and optional size adjustments
 *  used when computing molecular quadrature weights.
 */
struct MolecularWeightsSettings { 
    XCWeightAlg weight_alg = XCWeightAlg::SSF; ///< Weight partitioning scheme.
    bool becke_size_adjustment = false;        ///< Whether to use Becke size adjustments.
};


/**
 *  @brief Applies molecular partition weights to pre-generated quadrature tasks.
 *
 *  This class computes and modifies the quadrature weights stored in a
 *  LoadBalancer according to a specified partitioning scheme (e.g., SSF, Becke).
 */
class MolecularWeights {

public:

  using load_balancer_type = LoadBalancer;              ///< LoadBalancer type alias.
  using load_balancer_reference = load_balancer_type&;  ///< Reference to LoadBalancer.

private:

  using pimpl_type = detail::MolecularWeightsImpl;       ///< Implementation type.
  using pimpl_ptr_type = std::unique_ptr<pimpl_type>;    ///< Unique pointer to impl.
  pimpl_ptr_type pimpl_;                                 ///< Pointer to implementation instance.

public:

  /// Deleted default constructor.
  MolecularWeights() = delete;

  /// Destructor.
  ~MolecularWeights() noexcept;

  /**
   *  @brief Construct from a pre-built implementation.
   *  @param ptr Unique pointer to implementation instance.
   */
  MolecularWeights( pimpl_ptr_type&& ptr );

  /// Deleted copy constructor.
  MolecularWeights( const MolecularWeights& ) = delete;

  /// Move constructor.
  MolecularWeights( MolecularWeights&& ) noexcept;

  /**
   *  @brief Apply weight partitioning to local quadrature tasks.
   *  @param lb Reference to the LoadBalancer whose weights will be modified.
   */
  void modify_weights(load_balancer_reference lb) const;

  /**
   *  @brief Get the local timing tracker.
   *  @return Const reference to internal timer.
   */
  const util::Timer& get_timings() const;

}; // class MolecularWeights


/**
 *  @brief Factory for constructing MolecularWeights instances.
 *
 *  Provides a configurable interface to generate MolecularWeights objects
 *  for host or device execution spaces with specified settings.
 */
class MolecularWeightsFactory {

public:

    /// Deleted default constructor.
    MolecularWeightsFactory() = delete;

    /**
     *  @brief Construct a factory for generating MolecularWeights instances.
     *
     *  @param ex Execution space for weight computation.
     *            Acceptable values:
     *            - Host: Run on CPU.
     *            - Device: Run on GPU (if enabled).
     *  @param local_work_kernel_name Specification of the LocalWorkDriver kernel.
     *                                See LocalWorkDriver documentation for details.
     *  @param s Settings for the MolecularWeights calculation.
     */
    MolecularWeightsFactory( ExecutionSpace ex, 
                             std::string local_work_kernel_name,
                             MolecularWeightsSettings s);


    /**
     *  @brief Generate a shared-pointer MolecularWeights instance.
     *  @return Shared pointer to a new MolecularWeights object.
     */
    std::shared_ptr<MolecularWeights> get_shared_instance();

    /**
     *  @brief Generate a MolecularWeights instance.
     *  @return Newly constructed MolecularWeights object.
     */
    inline MolecularWeights get_instance(){
      return MolecularWeights( std::move( *get_shared_instance() ) );
    };

private:

    ExecutionSpace ex_;                  ///< Execution space for weight computation.
    std::string lwd_kernel_;             ///< LocalWorkDriver kernel name.
    MolecularWeightsSettings settings_;  ///< Settings for weight calculation.

}; // class MolecularWeightsFactory

}
