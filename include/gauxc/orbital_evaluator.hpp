/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy).
 *
 * (c) 2024-2026, Microsoft Corporation
 *
 * All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once

#include <cstdint>
#include <memory>

#include <gauxc/basisset.hpp>
#include <gauxc/enums.hpp>
#include <gauxc/external/cube.hpp>

namespace GauXC {

/** @brief Evaluate molecular orbitals and densities on arbitrary point sets.
 *
 *  Wraps the host collocation kernel exposed by the local work driver with a
 *  thread-parallel batched evaluation loop, hiding the driver factory and
 *  the per-thread AO scratch from callers. The class is designed to be
 *  constructed once per molecule and reused across many evaluations (e.g.
 *  one per active-space orbital) without re-initialising the driver.
 *
 *  The evaluator is intentionally decoupled from any particular file format:
 *  it returns plain numerical arrays. For Gaussian cube file I/O, see
 *  `gauxc/external/cube.hpp`.
 *
 *  Currently only `ExecutionSpace::Host` is supported; passing any other
 *  value will throw on construction. Device support can be added later by
 *  routing through a device-side collocation driver while keeping the same
 *  signatures.
 */
class OrbitalEvaluator {
 public:
  /** @brief Construct an evaluator bound to a basis set.
   *
   *  @param basis Basis set (copied internally).
   *  @param exec  Execution space; only `ExecutionSpace::Host` is supported.
   */
  explicit OrbitalEvaluator(BasisSet<double> basis,
                            ExecutionSpace exec = ExecutionSpace::Host);

  ~OrbitalEvaluator() noexcept;

  // Non-copyable, movable
  OrbitalEvaluator(const OrbitalEvaluator&) = delete;
  OrbitalEvaluator& operator=(const OrbitalEvaluator&) = delete;
  OrbitalEvaluator(OrbitalEvaluator&&) noexcept;
  OrbitalEvaluator& operator=(OrbitalEvaluator&&) noexcept;

  /// Number of basis functions (rows of the AO matrix).
  int32_t nbf() const noexcept;

  /// Underlying basis set.
  const BasisSet<double>& basis() const noexcept;

  /** @brief Evaluate a single MO chi(r) = sum_mu C[mu] * phi_mu(r).
   *
   *  @param[in]  npts    Number of evaluation points.
   *  @param[in]  points  AoS array of length 3*npts: (x0,y0,z0,x1,y1,z1,...)
   *                      in atomic units (Bohr).
   *  @param[in]  C       MO coefficient vector, length nbf().
   *  @param[out] out     Length-npts array of MO values.
   */
  void eval_orbital(int64_t npts, const double* points,
                    const double* C, double* out) const;

  /** @brief Evaluate `nmo` MOs simultaneously.
   *
   *  Storage layout:
   *    - `C`   : (nbf, nmo) column-major, leading dimension `ldc` (>= nbf).
   *    - `out` : (npts, nmo) column-major, leading dimension `ldo` (>= npts).
   *
   *  Equivalent to calling `eval_orbital` `nmo` times but amortises the AO
   *  collocation evaluation across all MOs (single AO buffer, GEMM contraction).
   */
  void eval_orbitals(int64_t npts, const double* points,
                     int32_t nmo, const double* C, int64_t ldc,
                     double* out, int64_t ldo) const;

  /** @brief Evaluate the electron density
   *         rho(r) = sum_{mu,nu} D[mu,nu] * phi_mu(r) * phi_nu(r).
   *
   *  @param[in]  npts    Number of evaluation points.
   *  @param[in]  points  AoS array of length 3*npts in Bohr.
   *  @param[in]  D       (nbf, nbf) symmetric density matrix, column-major,
   *                      leading dimension `ldd` (>= nbf).
   *  @param[out] out     Length-npts array of density values.
   */
  void eval_density(int64_t npts, const double* points,
                    const double* D, int64_t ldd,
                    double* out) const;

  /** @brief Evaluate a single MO on a CubeGrid without materialising all 3*N
   *         grid-point coordinates.
   */
  void eval_orbital(const CubeGrid& grid,
                    const double* C, double* out) const;

  /** @brief Evaluate `nmo` MOs on a CubeGrid without materialising all 3*N
   *         grid-point coordinates.
   */
  void eval_orbitals(const CubeGrid& grid,
                     int32_t nmo, const double* C, int64_t ldc,
                     double* out, int64_t ldo) const;

  /** @brief Evaluate the electron density on a CubeGrid without materialising
   *         all 3*N grid-point coordinates.
   */
  void eval_density(const CubeGrid& grid,
                    const double* D, int64_t ldd,
                    double* out) const;

 private:
  struct Impl;
  std::unique_ptr<Impl> pimpl_;
};

}  // namespace GauXC
