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

#include <cstddef>
#include <cstdint>
#include <memory>

#include <gauxc/basisset.hpp>
#include <gauxc/cube_grid.hpp>
#include <gauxc/enums.hpp>

namespace GauXC {

namespace detail {
  /// OrbitalEvaluator Implementation class
  class OrbitalEvaluatorImpl;
}

/** @brief Evaluate molecular orbitals and densities on arbitrary point sets.
 *
 *  Wraps the collocation kernel exposed by the local work driver with a
 *  thread-parallel batched evaluation loop, hiding the driver factory and
 *  the per-thread AO scratch from callers. The class is designed to be
 *  constructed once per molecule and reused across many evaluations (e.g.
 *  one per active-space orbital) without re-initialising the driver.
 *
 *  The evaluator is intentionally decoupled from any particular file format:
 *  it returns plain numerical arrays. For Gaussian cube file I/O, see
 *  `gauxc/external/cube.hpp`.
 *
 *  Points are evaluated in batches and each batch is screened against the
 *  per-shell cutoff radii. The batch decomposition is derived from the point
 *  count and the basis alone, never from the thread count, so results are
 *  bitwise reproducible from run to run and across thread counts.
 *
 *  The CubeGrid overloads are not bitwise equal to the pointer overloads given
 *  the same points, however. A CubeGrid is walked in spatially compact tiles,
 *  which screen better, whereas a caller-supplied point array is batched in
 *  the order it arrives, so the two drop different shells. Each neglected
 *  shell contributes at most the shell tolerance, which makes the difference
 *  of the order of the tolerance rather than bounded by it: it grows with the
 *  number of shells dropped, and reaches roughly twice the tolerance on
 *  benzene in cc-pVDZ.
 *
 *  Screening cost is governed by the shell tolerance, which the evaluator
 *  applies to its own private copy of the basis, so the caller's basis (and
 *  any SCF setup sharing it) is left untouched. Orbital error tracks the
 *  tolerance directly, density error goes as its square, so density tolerates
 *  a looser setting for the same accuracy: on a 110-atom system 1e-6 instead
 *  of the 1e-10 default runs ~2.5x faster on density.
 *
 *  Instances are produced by OrbitalEvaluatorFactory, which selects the
 *  implementation for a given ExecutionSpace.
 */
class OrbitalEvaluator {

  using pimpl_type = detail::OrbitalEvaluatorImpl;
  using pimpl_ptr_type = std::unique_ptr<pimpl_type>;
  pimpl_ptr_type pimpl_; ///< Pointer to implementation instance

public:

  // Delete default ctor
  OrbitalEvaluator() = delete;

  /// Construct an OrbitalEvaluator instance from a preconstructed implementation
  OrbitalEvaluator( pimpl_ptr_type&& pimpl );

  ~OrbitalEvaluator() noexcept;

  // Non-copyable, movable
  OrbitalEvaluator( const OrbitalEvaluator& ) = delete;
  OrbitalEvaluator& operator=( const OrbitalEvaluator& ) = delete;
  OrbitalEvaluator( OrbitalEvaluator&& ) noexcept;
  OrbitalEvaluator& operator=( OrbitalEvaluator&& ) noexcept;

  /// Number of basis functions (rows of the AO matrix).
  int32_t nbf() const;

  /// Underlying basis set.
  const BasisSet<double>& basis() const;

  /** @brief Evaluate a single MO chi(r) = sum_mu C[mu] * phi_mu(r).
   *
   *  @param[in]  npts    Number of evaluation points.
   *  @param[in]  points  AoS array of length 3*npts: (x0,y0,z0,x1,y1,z1,...)
   *                      in atomic units (Bohr).
   *  @param[in]  C       MO coefficient vector, length nbf().
   *  @param[out] out     Length-npts array of MO values.
   */
  void eval_orbital( size_t npts, const double* points,
                     const double* C, double* out ) const;

  /** @brief Evaluate `nmo` MOs simultaneously.
   *
   *  Storage layout:
   *    - `C`   : (nbf, nmo) column-major, leading dimension `ldc` (>= nbf).
   *    - `out` : (npts, nmo) column-major, leading dimension `ldo` (>= npts).
   *
   *  Equivalent to calling `eval_orbital` `nmo` times but amortises the AO
   *  collocation evaluation across all MOs (single AO buffer, GEMM contraction).
   */
  void eval_orbitals( size_t npts, const double* points,
                      int32_t nmo, const double* C, size_t ldc,
                      double* out, size_t ldo ) const;

  /** @brief Evaluate the electron density
   *         rho(r) = sum_{mu,nu} D[mu,nu] * phi_mu(r) * phi_nu(r).
   *
   *  @param[in]  npts    Number of evaluation points.
   *  @param[in]  points  AoS array of length 3*npts in Bohr.
   *  @param[in]  D       (nbf, nbf) symmetric density matrix, column-major,
   *                      leading dimension `ldd` (>= nbf).
   *  @param[out] out     Length-npts array of density values.
   */
  void eval_density( size_t npts, const double* points,
                     const double* D, size_t ldd,
                     double* out ) const;

  /** @brief Evaluate a single MO on a CubeGrid without materialising all 3*N
   *         grid-point coordinates.
   */
  void eval_orbital( const CubeGrid& grid,
                     const double* C, double* out ) const;

  /** @brief Evaluate `nmo` MOs on a CubeGrid without materialising all 3*N
   *         grid-point coordinates.
   */
  void eval_orbitals( const CubeGrid& grid,
                      int32_t nmo, const double* C, size_t ldc,
                      double* out, size_t ldo ) const;

  /** @brief Evaluate the electron density on a CubeGrid without materialising
   *         all 3*N grid-point coordinates.
   */
  void eval_density( const CubeGrid& grid,
                     const double* D, size_t ldd,
                     double* out ) const;

}; // class OrbitalEvaluator


/// A factory to generate OrbitalEvaluator instances
class OrbitalEvaluatorFactory {

public:

  // Delete default ctor
  OrbitalEvaluatorFactory() = delete;

  /**
   * @brief Construct an OrbitalEvaluator for a given execution space
   *
   * @param[in] ex    Execution space in which to evaluate orbitals/densities.
   *                  Currently only ExecutionSpace::Host is implemented;
   *                  anything else throws.
   * @param[in] basis Basis set (copied into the evaluator).
   * @param[in] screening_tolerance Shell tolerance applied to the evaluator's
   *                  own copy of `basis`, which sets the cutoff radii used
   *                  for per-batch screening. A basis carried in from an SCF
   *                  setup may be far tighter than cube-file precision needs
   *                  and costs several-fold here for no benefit.
   */
  static OrbitalEvaluator make_orbital_evaluator( ExecutionSpace ex,
    BasisSet<double> basis,
    double screening_tolerance = detail::default_shell_tolerance );

}; // class OrbitalEvaluatorFactory

}  // namespace GauXC
