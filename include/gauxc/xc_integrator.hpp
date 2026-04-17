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

#include <memory>

#include <gauxc/types.hpp>
#include <gauxc/load_balancer.hpp>
#include <gauxc/xc_integrator_settings.hpp>

namespace GauXC {

namespace detail {
  /// Forward declaration of XCIntegrator implementation.
  template <typename MatrixType>
  class XCIntegratorImpl;
}


/**
 *  @brief High-level interface for exchange-correlation and exact exchange integration.
 *
 *  XCIntegrator provides methods to evaluate XC energies, potentials, gradients,
 *  and exact exchange contributions using numerical quadrature on a molecular grid.
 *
 *  @tparam MatrixType Matrix type for density and Fock-like matrices.
 */
template <typename MatrixType>
class XCIntegrator {

public:

  using matrix_type   = MatrixType;                        ///< Matrix type alias.
  using value_type    = typename matrix_type::value_type;  ///< Scalar type alias.
  using basisset_type = BasisSet< value_type >;            ///< Basis set type alias.

  /// Return type for RKS exc+vxc: (energy, Vxc).
  using exc_vxc_type_rks  = std::tuple< value_type, matrix_type >;
  /// Return type for UKS exc+vxc: (energy, Vxc_alpha, Vxc_beta).
  using exc_vxc_type_uks  = std::tuple< value_type, matrix_type, matrix_type >;
  /// Return type for GKS exc+vxc: (energy, Vxc_scalar, Vxc_x, Vxc_y, Vxc_z).
  using exc_vxc_type_gks  = std::tuple< value_type, matrix_type, matrix_type, matrix_type, matrix_type >;
  /// Return type for XC gradient: vector of gradient components.
  using exc_grad_type = std::vector< value_type >;
  /// Return type for exact exchange matrix.
  using exx_type      = matrix_type;
  /// Return type for RKS Fxc contraction.
  using fxc_contraction_type_rks = matrix_type;
  /// Return type for UKS Fxc contraction: (Fxc_alpha, Fxc_beta).
  using fxc_contraction_type_uks = std::tuple< matrix_type, matrix_type >;
  /// Return type for dd-Psi evaluation.
  using dd_psi_type   = std::vector< value_type >;
  /// Return type for dd-Psi potential.
  using dd_psi_potential_type   = matrix_type;

private:

  using pimpl_type    = detail::XCIntegratorImpl<MatrixType>;  ///< Implementation type.

  std::unique_ptr<pimpl_type> pimpl_;  ///< PIMPL handle.

public:

  /// Default constructor.
  XCIntegrator() = default;

  /// Destructor.
  ~XCIntegrator() noexcept;

  /**
   *  @brief Construct from a pre-built implementation.
   *  @param pimpl Unique pointer to implementation instance.
   */
  XCIntegrator( std::unique_ptr<pimpl_type>&& pimpl );

  /// Deleted copy constructor.
  XCIntegrator( const XCIntegrator& ) = delete;

  /// Move constructor.
  XCIntegrator( XCIntegrator&& ) noexcept;

  /**
   *  @brief Integrate the electron density.
   *  @param P Density matrix.
   *  @return Integrated number of electrons.
   */
  value_type integrate_den( const MatrixType& P );

  /**
   *  @brief Evaluate the XC energy (RKS).
   *  @param P    Density matrix.
   *  @param opts Integration settings.
   *  @return XC energy.
   */
  value_type eval_exc( const MatrixType& P, const IntegratorSettingsXC& opts = IntegratorSettingsXC{} );

  /**
   *  @brief Evaluate the XC energy (UKS).
   *  @param Palpha Alpha density matrix.
   *  @param Pbeta  Beta density matrix.
   *  @param opts   Integration settings.
   *  @return XC energy.
   */
  value_type eval_exc( const MatrixType& Palpha, const MatrixType& Pbeta, const IntegratorSettingsXC& opts = IntegratorSettingsXC{} );

  /**
   *  @brief Evaluate the XC energy (GKS).
   *  @param Pscalar Scalar density matrix.
   *  @param Px      X-component spin density matrix.
   *  @param Py      Y-component spin density matrix.
   *  @param Pz      Z-component spin density matrix.
   *  @param opts    Integration settings.
   *  @return XC energy.
   */
  value_type eval_exc( const MatrixType& Pscalar, const MatrixType& Px, const MatrixType& Py, const MatrixType& Pz, const IntegratorSettingsXC& opts = IntegratorSettingsXC{} );

  /**
   *  @brief Evaluate XC energy and potential (RKS).
   *  @param P    Density matrix.
   *  @param opts Integration settings.
   *  @return Tuple of (energy, Vxc).
   */
  exc_vxc_type_rks eval_exc_vxc( const MatrixType& P, 
                                  const IntegratorSettingsXC& opts = IntegratorSettingsXC{} );

  /**
   *  @brief Evaluate XC energy and potential (UKS).
   *  @param Palpha Alpha density matrix.
   *  @param Pbeta  Beta density matrix.
   *  @param opts   Integration settings.
   *  @return Tuple of (energy, Vxc_alpha, Vxc_beta).
   */
  exc_vxc_type_uks eval_exc_vxc( const MatrixType& Palpha, const MatrixType& Pbeta,
                                  const IntegratorSettingsXC& opts = IntegratorSettingsXC{} );

  /**
   *  @brief Evaluate XC energy and potential (GKS).
   *  @param Pscalar Scalar density matrix.
   *  @param Px      X-component spin density matrix.
   *  @param Py      Y-component spin density matrix.
   *  @param Pz      Z-component spin density matrix.
   *  @param opts    Integration settings.
   *  @return Tuple of (energy, Vxc_scalar, Vxc_x, Vxc_y, Vxc_z).
   */
  exc_vxc_type_gks eval_exc_vxc( const MatrixType& Pscalar, const MatrixType& Px, const MatrixType& Py, const MatrixType& Pz,
                                  const IntegratorSettingsXC& opts = IntegratorSettingsXC{} );

  /**
   *  @brief Evaluate XC gradient (RKS).
   *  @param P    Density matrix.
   *  @param opts Integration settings.
   *  @return Vector of gradient components (3*natoms).
   */
  exc_grad_type eval_exc_grad( const MatrixType& P, const IntegratorSettingsXC& opts = IntegratorSettingsXC{} );

  /**
   *  @brief Evaluate XC gradient (UKS).
   *  @param Palpha Alpha density matrix.
   *  @param Pbeta  Beta density matrix.
   *  @param opts   Integration settings.
   *  @return Vector of gradient components (3*natoms).
   */
  exc_grad_type eval_exc_grad( const MatrixType& Palpha, const MatrixType& Pbeta, const IntegratorSettingsXC& opts = IntegratorSettingsXC{} );

  /**
   *  @brief Evaluate exact exchange matrix (sn-LinK).
   *  @param P    Density matrix.
   *  @param opts EXX integration settings.
   *  @return Exact exchange matrix K.
   */
  exx_type eval_exx( const MatrixType& P, 
                     const IntegratorSettingsEXX& opts = IntegratorSettingsEXX{} );

  /**
   *  @brief Evaluate Fxc contraction (RKS).
   *  @param P    Density matrix.
   *  @param X    Trial vector for contraction.
   *  @param opts Integration settings.
   *  @return Contracted Fxc matrix.
   */
  fxc_contraction_type_rks eval_fxc_contraction( const MatrixType& P, const MatrixType& X,
                                const IntegratorSettingsXC& opts = IntegratorSettingsXC{} );

  /**
   *  @brief Evaluate Fxc contraction (UKS).
   *  @param Palpha Alpha density matrix.
   *  @param Pbeta  Beta density matrix.
   *  @param Xalpha Alpha trial vector.
   *  @param Xbeta  Beta trial vector.
   *  @param opts   Integration settings.
   *  @return Tuple of (Fxc_alpha, Fxc_beta).
   */
  fxc_contraction_type_uks eval_fxc_contraction( const MatrixType& Palpha, const MatrixType& Pbeta, const MatrixType& Xalpha, const MatrixType& Xbeta,
                                const IntegratorSettingsXC& opts = IntegratorSettingsXC{} );

  /**
   *  @brief Evaluate dd-Psi diagnostic.
   *  @param P     Density matrix.
   *  @param deriv Derivative order.
   *  @return Vector of dd-Psi values.
   */
  dd_psi_type eval_dd_psi( const MatrixType& P, unsigned deriv );

  /**
   *  @brief Evaluate dd-Psi potential.
   *  @param P     Density matrix.
   *  @param deriv Derivative order.
   *  @return dd-Psi potential matrix.
   */
  dd_psi_potential_type eval_dd_psi_potential( const MatrixType& P, unsigned deriv );

  /**
   *  @brief Get the internal timing tracker.
   *  @return Const reference to timer.
   */
  const util::Timer& get_timings() const;

  /**
   *  @brief Get the associated load balancer (const).
   *  @return Const reference to LoadBalancer.
   */
  const LoadBalancer& load_balancer() const;

  /**
   *  @brief Get the associated load balancer (mutable).
   *  @return Reference to LoadBalancer.
   */
  LoadBalancer& load_balancer();
};


}
