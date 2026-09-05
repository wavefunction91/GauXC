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
#include "ut_common.hpp"
#include "catch2/catch.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <gauxc/external/cube.hpp>
#include <gauxc/orbital_evaluator.hpp>
#include <gauxc/xc_integrator/local_work_driver.hpp>
#include <gauxc/gauxc_config.hpp>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef GAUXC_HAS_HDF5
#include <highfive/H5File.hpp>
#endif
#ifdef GAUXC_HAS_MPI
#include <mpi.h>
#endif

#include "standards.hpp"

// Reference implementation lives in src/, behind the public header surface.
#include "xc_integrator/local_work_driver/host/local_host_work_driver.hpp"

using namespace GauXC;

namespace {

/// Unique path for a test output file, inside the build-tree scratch
/// directory configured by CMake.
std::string make_temp_path(const char* suffix) {
  static int counter = 0;
  return std::string(GAUXC_TEST_TMP_PATH) + "/gauxc_test_" +
         std::to_string(++counter) + suffix;
}

OrbitalEvaluator make_evaluator(const BasisSet<double>& basis, double tol) {
  return OrbitalEvaluatorFactory::make_orbital_evaluator(ExecutionSpace::Host,
                                                         basis, tol);
}

/// AO collocation over every shell, i.e. the unscreened reference.
std::vector<double> reference_collocation(const BasisSet<double>& basis,
                                          int64_t npts, const double* pts) {
  const int32_t nbf = basis.nbf();
  std::vector<double> ao(static_cast<size_t>(nbf) * npts);
  auto drv = LocalWorkDriverFactory::make_local_work_driver(
      ExecutionSpace::Host, "Reference");
  auto* host_drv = dynamic_cast<LocalHostWorkDriver*>(drv.get());
  REQUIRE(host_drv != nullptr);
  std::vector<int32_t> shell_list(basis.size());
  for (size_t i = 0; i < shell_list.size(); ++i)
    shell_list[i] = static_cast<int32_t>(i);
  host_drv->eval_collocation(static_cast<size_t>(npts),
                             static_cast<size_t>(basis.nshells()),
                             static_cast<size_t>(nbf), pts, basis,
                             shell_list.data(), ao.data());
  return ao;
}

/// Build a deterministic PRNG-based set of points within a small bounding box
/// around the molecule. Avoids putting samples too close to nuclei to keep
/// values numerically well-behaved.
std::vector<double> make_random_points(int64_t npts, unsigned seed = 1234u) {
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> u(-2.5, 2.5);
  std::vector<double> pts(static_cast<size_t>(npts) * 3);
  for (int64_t p = 0; p < npts; ++p) {
    pts[3 * p + 0] = u(gen);
    pts[3 * p + 1] = u(gen);
    pts[3 * p + 2] = u(gen);
  }
  return pts;
}

std::vector<double> make_random_vector(size_t n, unsigned seed) {
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> u(-1.0, 1.0);
  std::vector<double> v(n);
  for (auto& x : v) x = u(gen);
  return v;
}

}  // namespace

TEST_CASE("OrbitalEvaluator / Water cc-pVDZ matches eval_collocation",
          "[orbital_evaluator]") {
  auto mol = make_water();
  auto basis = make_ccpvdz(mol, SphericalType(true));
  for (auto& sh : basis) sh.set_shell_tolerance(1e-12);

  const int32_t nbf = basis.nbf();
  REQUIRE(nbf > 0);

  const int64_t npts = 137;  // not a multiple of any batch size
  const auto pts = make_random_points(npts);

  // Reference: AO collocation directly via the LocalHostWorkDriver.
  const auto ao_ref = reference_collocation(basis, npts, pts.data());

  auto eval = make_evaluator(basis, 1e-12);
  REQUIRE(eval.nbf() == nbf);

  SECTION("eval_orbital with one-hot coefficient reproduces single AO column") {
    std::vector<double> C(nbf, 0.0);
    std::vector<double> out(static_cast<size_t>(npts), 0.0);
    for (int32_t mu : {0, nbf / 3, nbf / 2, nbf - 1}) {
      std::fill(C.begin(), C.end(), 0.0);
      C[static_cast<size_t>(mu)] = 1.0;
      std::fill(out.begin(), out.end(), 0.0);
      eval.eval_orbital(npts, pts.data(), C.data(), out.data());
      for (int64_t p = 0; p < npts; ++p) {
        const double ref =
            ao_ref[static_cast<size_t>(p) * nbf + static_cast<size_t>(mu)];
        CHECK(out[static_cast<size_t>(p)] == Approx(ref).margin(1e-12));
      }
    }
  }

  SECTION("eval_orbitals with random C matches AO ^T @ C") {
    const int32_t nmo = 4;
    const auto C = make_random_vector(static_cast<size_t>(nbf) * nmo, 99u);

    std::vector<double> out(static_cast<size_t>(npts) * nmo, 0.0);
    eval.eval_orbitals(npts, pts.data(), nmo, C.data(), nbf, out.data(), npts);

    for (int32_t j = 0; j < nmo; ++j) {
      for (int64_t p = 0; p < npts; ++p) {
        double ref = 0.0;
        for (int32_t mu = 0; mu < nbf; ++mu) {
          ref += C[static_cast<size_t>(j) * nbf + mu] *
                 ao_ref[static_cast<size_t>(p) * nbf + mu];
        }
        const double got =
            out[static_cast<size_t>(j) * npts + static_cast<size_t>(p)];
        CHECK(got == Approx(ref).margin(1e-10));
      }
    }
  }

  SECTION("eval_density with identity D equals sum of squared AO values") {
    std::vector<double> D(static_cast<size_t>(nbf) * nbf, 0.0);
    for (int32_t mu = 0; mu < nbf; ++mu) {
      D[static_cast<size_t>(mu) * nbf + mu] = 1.0;
    }
    std::vector<double> out(static_cast<size_t>(npts), 0.0);
    eval.eval_density(npts, pts.data(), D.data(), nbf, out.data());

    for (int64_t p = 0; p < npts; ++p) {
      double ref = 0.0;
      for (int32_t mu = 0; mu < nbf; ++mu) {
        const double a = ao_ref[static_cast<size_t>(p) * nbf + mu];
        ref += a * a;
      }
      CHECK(out[static_cast<size_t>(p)] == Approx(ref).margin(1e-10));
    }
  }

  SECTION("eval_density with rank-1 D = c c^T equals (c.AO)^2") {
    const auto c = make_random_vector(static_cast<size_t>(nbf), 7u);

    std::vector<double> D(static_cast<size_t>(nbf) * nbf, 0.0);
    for (int32_t mu = 0; mu < nbf; ++mu) {
      for (int32_t nu = 0; nu < nbf; ++nu) {
        D[static_cast<size_t>(mu) * nbf + nu] = c[mu] * c[nu];
      }
    }
    std::vector<double> out(static_cast<size_t>(npts), 0.0);
    eval.eval_density(npts, pts.data(), D.data(), nbf, out.data());

    std::vector<double> orb(static_cast<size_t>(npts), 0.0);
    eval.eval_orbital(npts, pts.data(), c.data(), orb.data());

    for (int64_t p = 0; p < npts; ++p) {
      const double ref = orb[static_cast<size_t>(p)] * orb[static_cast<size_t>(p)];
      CHECK(out[static_cast<size_t>(p)] == Approx(ref).margin(1e-10));
    }
  }

  SECTION("invalid leading dimensions throw") {
    std::vector<double> C(nbf, 0.0);
    std::vector<double> out(static_cast<size_t>(npts), 0.0);
    CHECK_THROWS(eval.eval_orbitals(npts, pts.data(), 1, C.data(), nbf - 1,
                                    out.data(), npts));
    CHECK_THROWS(eval.eval_orbitals(npts, pts.data(), 1, C.data(), nbf,
                                    out.data(), npts - 1));
    // ldo is narrowed to int for the BLAS call; oversized values must throw
    // rather than silently truncate.
    const size_t huge_ldo =
        static_cast<size_t>(std::numeric_limits<int>::max()) + 1;
    CHECK_THROWS(eval.eval_orbitals(npts, pts.data(), 1, C.data(), nbf,
                                    out.data(), huge_ldo));
  }
}

TEST_CASE("OrbitalEvaluator public API surface", "[orbital_evaluator]") {
  auto mol = make_water();
  constexpr double shell_tol = 1e-10;
  auto basis = make_ccpvdz(mol, SphericalType(true));
  for (auto& sh : basis) sh.set_shell_tolerance(shell_tol);
  const int32_t nbf = basis.nbf();

  SECTION("basis() exposes the evaluator's own retuned copy") {
    auto tight = make_ccpvdz(mol, SphericalType(true));
    for (auto& sh : tight) sh.set_shell_tolerance(1e-12);
    auto eval = make_evaluator(tight, 1e-4);

    REQUIRE(eval.basis().nbf() == nbf);
    REQUIRE(eval.basis().size() == tight.size());

    // The caller's basis keeps the tolerance it was given, so the evaluator's
    // looser one must show up as shorter cutoff radii on its copy alone.
    bool any_shorter = false;
    for (size_t i = 0; i < tight.size(); ++i) {
      CHECK(eval.basis()[i].cutoff_radius() <= tight[i].cutoff_radius());
      if (eval.basis()[i].cutoff_radius() < tight[i].cutoff_radius())
        any_shorter = true;
    }
    CHECK(any_shorter);
  }

  SECTION("move construction and assignment carry the implementation") {
    const int64_t npts = 64;
    const auto pts = make_random_points(npts, 7u);
    const auto C = make_random_vector(static_cast<size_t>(nbf), 11u);

    auto eval = make_evaluator(basis, shell_tol);
    std::vector<double> ref(static_cast<size_t>(npts));
    eval.eval_orbital(npts, pts.data(), C.data(), ref.data());

    OrbitalEvaluator moved(std::move(eval));
    REQUIRE(moved.nbf() == nbf);
    std::vector<double> out(static_cast<size_t>(npts));
    moved.eval_orbital(npts, pts.data(), C.data(), out.data());
    for (int64_t p = 0; p < npts; ++p) CHECK(out[p] == ref[p]);

    // Assigned over an evaluator built with a looser tolerance, so the result
    // would shift if the target's original implementation survived.
    auto target = make_evaluator(basis, 1e-3);
    target = std::move(moved);
    REQUIRE(target.nbf() == nbf);
    std::vector<double> out2(static_cast<size_t>(npts));
    target.eval_orbital(npts, pts.data(), C.data(), out2.data());
    for (int64_t p = 0; p < npts; ++p) CHECK(out2[p] == ref[p]);
  }

  SECTION("a non-Host execution space is rejected") {
    CHECK_THROWS(OrbitalEvaluatorFactory::make_orbital_evaluator(
        ExecutionSpace::Device, basis));
  }

  SECTION("null pointers and undersized leading dimensions throw") {
    auto eval = make_evaluator(basis, shell_tol);
    const int64_t npts = 8;
    const auto pts = make_random_points(npts, 3u);
    std::vector<double> C(static_cast<size_t>(nbf), 1.0);
    std::vector<double> out(static_cast<size_t>(npts), 0.0);
    std::vector<double> D(static_cast<size_t>(nbf) * nbf, 0.0);

    CHECK_THROWS(eval.eval_orbital(npts, nullptr, C.data(), out.data()));
    CHECK_THROWS(eval.eval_orbital(npts, pts.data(), nullptr, out.data()));
    CHECK_THROWS(eval.eval_orbital(npts, pts.data(), C.data(), nullptr));
    CHECK_THROWS(eval.eval_density(npts, nullptr, D.data(), nbf, out.data()));
    CHECK_THROWS(eval.eval_density(npts, pts.data(), nullptr, nbf, out.data()));
    CHECK_THROWS(eval.eval_density(npts, pts.data(), D.data(), nbf, nullptr));
    CHECK_THROWS(
        eval.eval_density(npts, pts.data(), D.data(), nbf - 1, out.data()));

    const auto grid = CubeGrid::from_molecule(mol, 3, 3, 3);
    std::vector<double> gout(static_cast<size_t>(grid.num_points()), 0.0);
    CHECK_THROWS(eval.eval_orbital(grid, nullptr, gout.data()));
    CHECK_THROWS(eval.eval_density(grid, D.data(), nbf - 1, gout.data()));
  }

  SECTION("empty work is a no-op rather than an error") {
    auto eval = make_evaluator(basis, shell_tol);
    std::vector<double> C(static_cast<size_t>(nbf), 1.0);
    std::vector<double> D(static_cast<size_t>(nbf) * nbf, 0.0);
    constexpr double untouched = -7.0;
    std::vector<double> out(4, untouched);
    const auto pts = make_random_points(4, 5u);

    CHECK_NOTHROW(
        eval.eval_orbitals(0, nullptr, 1, C.data(), nbf, out.data(), 0));
    CHECK_NOTHROW(eval.eval_density(0, nullptr, D.data(), nbf, out.data()));
    CHECK_NOTHROW(
        eval.eval_orbitals(4, pts.data(), 0, C.data(), nbf, out.data(), 4));

    CubeGrid empty;
    empty.nx = 0;
    CHECK_NOTHROW(eval.eval_orbital(empty, C.data(), out.data()));
    CHECK_NOTHROW(eval.eval_density(empty, D.data(), nbf, out.data()));

    for (double v : out) CHECK(v == untouched);
  }
}

TEST_CASE("CubeGrid eval overloads match pointer-based eval",
          "[orbital_evaluator]") {
  auto mol = make_water();
  auto basis = make_ccpvdz(mol, SphericalType(true));
  for (auto& sh : basis) sh.set_shell_tolerance(1e-12);
  const int32_t nbf = basis.nbf();
  auto eval = make_evaluator(basis, 1e-12);

  // 6x7x8 fits in a single batch; 20x20x20 = 8000 points spreads over many
  // batches at any thread count, exercising the trailing partial batch and
  // the dynamic schedule.
  const std::vector<CubeGrid> grids = {CubeGrid::from_molecule(mol, 6, 7, 8),
                                       CubeGrid::from_molecule(mol, 20, 20, 20)};

  const auto C = make_random_vector(static_cast<size_t>(nbf), 42u);

  SECTION("eval_orbital(grid) matches eval_orbital(npts, points)") {
    for (const auto& grid : grids) {
      const int64_t npts = grid.num_points();
      const auto pts = grid.points();

      std::vector<double> ref(static_cast<size_t>(npts));
      eval.eval_orbital(npts, pts.data(), C.data(), ref.data());

      std::vector<double> out(static_cast<size_t>(npts));
      eval.eval_orbital(grid, C.data(), out.data());

      for (int64_t p = 0; p < npts; ++p) {
        CHECK(out[p] == Approx(ref[p]).margin(1e-12));
      }
    }
  }

  SECTION("eval_density(grid) matches eval_density(npts, points)") {
    // Identity density.
    std::vector<double> D(static_cast<size_t>(nbf) * nbf, 0.0);
    for (int32_t i = 0; i < nbf; ++i) D[i * nbf + i] = 1.0;

    for (const auto& grid : grids) {
      const int64_t npts = grid.num_points();
      const auto pts = grid.points();

      std::vector<double> ref(static_cast<size_t>(npts));
      eval.eval_density(npts, pts.data(), D.data(), nbf, ref.data());

      std::vector<double> out(static_cast<size_t>(npts));
      eval.eval_density(grid, D.data(), nbf, out.data());

      for (int64_t p = 0; p < npts; ++p) {
        CHECK(out[p] == Approx(ref[p]).margin(1e-12));
      }
    }
  }
}

TEST_CASE("OrbitalEvaluator shell screening", "[orbital_evaluator]") {
  // Two well-separated centres of different elements: batches near one centre
  // screen out every shell of the other, so 0 < nbe < nbf is reached, and
  // batches in the gap screen out everything (nbe == 0).
  Molecule mol;
  mol.emplace_back(AtomicNumber(8), 0.0, 0.0, 0.0);
  mol.emplace_back(AtomicNumber(1), 20.0, 0.0, 0.0);

  auto basis = make_ccpvdz(mol, SphericalType(true));
  constexpr double shell_tol = 1e-10;
  for (auto& sh : basis) sh.set_shell_tolerance(shell_tol);
  const int32_t nbf = basis.nbf();
  auto eval = make_evaluator(basis, shell_tol);

  const auto C = make_random_vector(static_cast<size_t>(nbf), 2024u);
  std::vector<double> D(static_cast<size_t>(nbf) * nbf, 0.0);
  for (int32_t i = 0; i < nbf; ++i) D[i * nbf + i] = 1.0;

  SECTION("a grid far from every atom evaluates to exactly zero") {
    CubeGrid far_grid;
    far_grid.origin = {100.0, 100.0, 100.0};
    far_grid.spacing = {0.5, 0.5, 0.5};
    far_grid.nx = far_grid.ny = far_grid.nz = 12;
    const int64_t npts = far_grid.num_points();

    std::vector<double> orb(static_cast<size_t>(npts), 1.0);
    eval.eval_orbital(far_grid, C.data(), orb.data());
    for (int64_t p = 0; p < npts; ++p) CHECK(orb[p] == 0.0);

    std::vector<double> rho(static_cast<size_t>(npts), 1.0);
    eval.eval_density(far_grid, D.data(), nbf, rho.data());
    for (int64_t p = 0; p < npts; ++p) CHECK(rho[p] == 0.0);
  }

  SECTION("partial screening agrees with the unscreened reference") {
    // Grid elongated along x (the slow axis), so a batch of consecutive
    // points is a narrow x-slab: near one centre, near the other, or in the
    // empty gap between them.
    CubeGrid grid;
    grid.origin = {-4.0, -4.0, -4.0};
    grid.spacing = {0.8, 2.0, 2.0};
    grid.nx = 40;
    grid.ny = 4;
    grid.nz = 4;
    const int64_t npts = grid.num_points();
    const auto pts = grid.points();
    const auto ao_ref = reference_collocation(basis, npts, pts.data());

    std::vector<double> orb(static_cast<size_t>(npts));
    eval.eval_orbital(grid, C.data(), orb.data());

    std::vector<double> rho(static_cast<size_t>(npts));
    eval.eval_density(grid, D.data(), nbf, rho.data());

    // The grid overload walks spatial tiles while the pointer overload takes
    // contiguous index ranges, so the two screen against different bounding
    // boxes. They agree to within the shell tolerance, not bitwise.
    std::vector<double> orb_pts(static_cast<size_t>(npts));
    eval.eval_orbital(npts, pts.data(), C.data(), orb_pts.data());
    std::vector<double> rho_pts(static_cast<size_t>(npts));
    eval.eval_density(npts, pts.data(), D.data(), nbf, rho_pts.data());

    bool any_nonzero = false;
    for (int64_t p = 0; p < npts; ++p) {
      CHECK(orb[p] == Approx(orb_pts[p]).margin(shell_tol));
      CHECK(rho[p] == Approx(rho_pts[p]).margin(shell_tol));

      double orb_ref = 0.0, rho_ref = 0.0;
      for (int32_t mu = 0; mu < nbf; ++mu) {
        const double a = ao_ref[static_cast<size_t>(p) * nbf + mu];
        orb_ref += C[static_cast<size_t>(mu)] * a;
        rho_ref += a * a;
      }
      if (std::fabs(orb_ref) > 1e-3) any_nonzero = true;
      // Screening discards shells whose magnitude is below the shell
      // tolerance across the batch bbox, so agreement is at that level.
      CHECK(orb[p] == Approx(orb_ref).margin(1e-9));
      CHECK(rho[p] == Approx(rho_ref).margin(1e-9));
    }
    CHECK(any_nonzero);
  }
}

TEST_CASE("OrbitalEvaluator tiled grid traversal matches the reference",
          "[orbital_evaluator]") {
  // Grids are walked in spatial tiles rather than contiguous index ranges, so
  // a tile's points are scattered into the output. The two centres here are
  // far enough apart that screening is active, and the 32 Bohr z extent is
  // several times the largest cc-pVDZ cutoff radius (13.2 Bohr), so tiles at
  // opposite ends of the grid see different shell sets.
  Molecule mol;
  mol.emplace_back(AtomicNumber(8), 0.0, 0.0, 0.0);
  mol.emplace_back(AtomicNumber(1), 0.0, 0.0, 20.0);

  constexpr double shell_tol = 1e-10;
  auto basis = make_ccpvdz(mol, SphericalType(true));
  for (auto& sh : basis) sh.set_shell_tolerance(shell_tol);
  const int32_t nbf = basis.nbf();
  auto eval = make_evaluator(basis, shell_tol);

  CubeGrid grid;
  grid.origin = {-4.0, -4.0, -6.0};
  grid.spacing = {2.0, 2.0, 0.8};
  grid.nx = 4;
  grid.ny = 4;
  grid.nz = 40;
  const int64_t npts = grid.num_points();
  const auto pts = grid.points();
  const auto ao_ref = reference_collocation(basis, npts, pts.data());

  const auto C = make_random_vector(static_cast<size_t>(nbf), 8675u);
  std::vector<double> D(static_cast<size_t>(nbf) * nbf, 0.0);
  for (int32_t i = 0; i < nbf; ++i) D[i * nbf + i] = 1.0;

  std::vector<double> orb(static_cast<size_t>(npts));
  std::vector<double> rho(static_cast<size_t>(npts));
  eval.eval_orbital(grid, C.data(), orb.data());
  eval.eval_density(grid, D.data(), nbf, rho.data());

  // The pointer overload always uses contiguous batches, so it cross-checks
  // the tiled result against a different decomposition.
  std::vector<double> orb_pts(static_cast<size_t>(npts));
  std::vector<double> rho_pts(static_cast<size_t>(npts));
  eval.eval_orbital(npts, pts.data(), C.data(), orb_pts.data());
  eval.eval_density(npts, pts.data(), D.data(), nbf, rho_pts.data());

  bool any_nonzero = false;
  for (int64_t p = 0; p < npts; ++p) {
    double orb_ref = 0.0, rho_ref = 0.0;
    for (int32_t mu = 0; mu < nbf; ++mu) {
      const double a = ao_ref[static_cast<size_t>(p) * nbf + mu];
      orb_ref += C[static_cast<size_t>(mu)] * a;
      rho_ref += a * a;
    }
    if (std::fabs(orb_ref) > 1e-3) any_nonzero = true;
    CHECK(orb[p] == Approx(orb_ref).margin(1e-9));
    CHECK(rho[p] == Approx(rho_ref).margin(1e-9));
    CHECK(orb[p] == Approx(orb_pts[p]).margin(shell_tol));
    CHECK(rho[p] == Approx(rho_pts[p]).margin(shell_tol));
  }
  CHECK(any_nonzero);
}

TEST_CASE("OrbitalEvaluator scatters multiple orbitals into a padded output",
          "[orbital_evaluator]") {
  // A tiled batch is not contiguous in the output, so its result is staged and
  // scattered column by column. nmo > 1 with a padded ldo is the only
  // combination that exercises both strides of that scatter. Padding entries
  // hold sentinels: C must never be read past nbf, out never written past npts.
  Molecule mol;
  mol.emplace_back(AtomicNumber(8), 0.0, 0.0, 0.0);
  mol.emplace_back(AtomicNumber(1), 0.0, 0.0, 20.0);

  constexpr double shell_tol = 1e-10;
  auto basis = make_ccpvdz(mol, SphericalType(true));
  for (auto& sh : basis) sh.set_shell_tolerance(shell_tol);
  const int32_t nbf = basis.nbf();
  auto eval = make_evaluator(basis, shell_tol);

  CubeGrid grid;
  grid.origin = {-4.0, -4.0, -6.0};
  grid.spacing = {2.0, 2.0, 0.8};
  grid.nx = 4;
  grid.ny = 4;
  grid.nz = 40;
  const int64_t npts = grid.num_points();
  const auto pts = grid.points();
  const auto ao_ref = reference_collocation(basis, npts, pts.data());

  constexpr int32_t nmo = 3;
  constexpr double sentinel = -12345.0;
  const size_t ldc = static_cast<size_t>(nbf) + 7;
  const size_t ldo = static_cast<size_t>(npts) + 5;

  const auto Craw = make_random_vector(static_cast<size_t>(nbf) * nmo, 99u);
  std::vector<double> C(ldc * nmo, sentinel);
  for (int32_t j = 0; j < nmo; ++j)
    for (int32_t mu = 0; mu < nbf; ++mu)
      C[j * ldc + mu] = Craw[static_cast<size_t>(j) * nbf + mu];

  std::vector<double> out(ldo * nmo, sentinel);
  eval.eval_orbitals(grid, nmo, C.data(), ldc, out.data(), ldo);

  // The pointer overload always batches contiguously, so it checks the tiled
  // scatter against a different decomposition of the same grid.
  std::vector<double> out_pts(ldo * nmo, sentinel);
  eval.eval_orbitals(npts, pts.data(), nmo, C.data(), ldc, out_pts.data(), ldo);

  for (int32_t j = 0; j < nmo; ++j) {
    for (int64_t p = 0; p < npts; ++p) {
      double ref = 0.0;
      for (int32_t mu = 0; mu < nbf; ++mu)
        ref += Craw[static_cast<size_t>(j) * nbf + mu] *
               ao_ref[static_cast<size_t>(p) * nbf + mu];
      const size_t k = static_cast<size_t>(j) * ldo + static_cast<size_t>(p);
      CHECK(out[k] == Approx(ref).margin(1e-9));
      CHECK(out[k] == Approx(out_pts[k]).margin(shell_tol));
    }
    for (size_t p = static_cast<size_t>(npts); p < ldo; ++p) {
      CHECK(out[static_cast<size_t>(j) * ldo + p] == sentinel);
      CHECK(out_pts[static_cast<size_t>(j) * ldo + p] == sentinel);
    }
  }
}

TEST_CASE("OrbitalEvaluator handles single-plane grids", "[orbital_evaluator]") {
  // An axis of one point gets zero spacing from from_molecule, which the
  // tile-shaping code must treat as unconstrained rather than divide by.
  auto mol = make_water();
  constexpr double shell_tol = 1e-10;
  auto basis = make_ccpvdz(mol, SphericalType(true));
  for (auto& sh : basis) sh.set_shell_tolerance(shell_tol);
  const int32_t nbf = basis.nbf();
  auto eval = make_evaluator(basis, shell_tol);

  const auto C = make_random_vector(static_cast<size_t>(nbf), 4242u);
  std::vector<double> D(static_cast<size_t>(nbf) * nbf, 0.0);
  for (int32_t i = 0; i < nbf; ++i) D[static_cast<size_t>(i) * nbf + i] = 1.0;

  const std::vector<CubeGrid> grids = {CubeGrid::from_molecule(mol, 1, 9, 9),
                                       CubeGrid::from_molecule(mol, 9, 1, 9),
                                       CubeGrid::from_molecule(mol, 9, 9, 1)};

  for (const auto& grid : grids) {
    const int64_t npts = grid.num_points();
    const auto pts = grid.points();
    const auto ao = reference_collocation(basis, npts, pts.data());

    std::vector<double> orb(static_cast<size_t>(npts));
    std::vector<double> rho(static_cast<size_t>(npts));
    eval.eval_orbital(grid, C.data(), orb.data());
    eval.eval_density(grid, D.data(), nbf, rho.data());

    for (int64_t p = 0; p < npts; ++p) {
      double orb_ref = 0.0, rho_ref = 0.0;
      for (int32_t mu = 0; mu < nbf; ++mu) {
        const double a = ao[static_cast<size_t>(p) * nbf + mu];
        orb_ref += C[static_cast<size_t>(mu)] * a;
        rho_ref += a * a;
      }
      CHECK(orb[p] == Approx(orb_ref).margin(1e-9));
      CHECK(rho[p] == Approx(rho_ref).margin(1e-9));
    }
  }
}

TEST_CASE("OrbitalEvaluator density with a general D and padded ldd",
          "[orbital_evaluator]") {
  // Every other density test uses a diagonal D and ldd == nbf, and neither
  // can see a whole class of defect. Permuting an identity on both sides
  // leaves it unchanged, so a scrambled compressed-submatrix map is invisible
  // to it; and ldd == nbf makes a leading dimension dropped in favour of nbf
  // a no-op. Screening has to be active for the compressed submatrix to
  // differ from D at all, hence the two well-separated centres.
  Molecule mol;
  mol.emplace_back(AtomicNumber(8), 0.0, 0.0, 0.0);
  mol.emplace_back(AtomicNumber(8), 0.0, 0.0, 20.0);

  constexpr double shell_tol = 1e-10;
  auto basis = make_ccpvdz(mol, SphericalType(true));
  for (auto& sh : basis) sh.set_shell_tolerance(shell_tol);
  const int32_t nbf = basis.nbf();
  auto eval = make_evaluator(basis, shell_tol);

  CubeGrid grid;
  grid.origin = {-3.0, -3.0, -4.0};
  grid.spacing = {2.0, 2.0, 1.0};
  grid.nx = 4;
  grid.ny = 4;
  grid.nz = 28;
  const int64_t npts = grid.num_points();
  const auto pts = grid.points();
  const auto ao = reference_collocation(basis, npts, pts.data());

  constexpr double sentinel = -4321.0;
  const size_t ldd = static_cast<size_t>(nbf) + 6;

  const auto rnd = make_random_vector(static_cast<size_t>(nbf) * nbf, 31337u);
  std::vector<double> Dsq(static_cast<size_t>(nbf) * nbf);
  for (int32_t i = 0; i < nbf; ++i)
    for (int32_t j = 0; j < nbf; ++j)
      Dsq[static_cast<size_t>(j) * nbf + i] =
          rnd[static_cast<size_t>(j) * nbf + i] +
          rnd[static_cast<size_t>(i) * nbf + j];

  // Rows past nbf are poisoned, so an ldd mistaken for nbf reads them.
  std::vector<double> D(ldd * nbf, sentinel);
  for (int32_t j = 0; j < nbf; ++j)
    for (int32_t i = 0; i < nbf; ++i)
      D[static_cast<size_t>(j) * ldd + i] = Dsq[static_cast<size_t>(j) * nbf + i];

  std::vector<double> rho_grid(static_cast<size_t>(npts));
  std::vector<double> rho_pts(static_cast<size_t>(npts));
  eval.eval_density(grid, D.data(), ldd, rho_grid.data());
  eval.eval_density(npts, pts.data(), D.data(), ldd, rho_pts.data());

  double max_ref = 0.0;
  for (int64_t p = 0; p < npts; ++p) {
    const double* a = ao.data() + static_cast<size_t>(p) * nbf;
    double ref = 0.0;
    for (int32_t mu = 0; mu < nbf; ++mu) {
      double acc = 0.0;
      for (int32_t nu = 0; nu < nbf; ++nu)
        acc += Dsq[static_cast<size_t>(nu) * nbf + mu] * a[nu];
      ref += acc * a[mu];
    }
    max_ref = std::max(max_ref, std::fabs(ref));
    CHECK(rho_grid[p] == Approx(ref).margin(1e-8));
    CHECK(rho_pts[p] == Approx(ref).margin(1e-8));
  }
  CHECK(max_ref > 1e-2);
}

TEST_CASE("OrbitalEvaluator boxes a batch out to its far face",
          "[orbital_evaluator]") {
  // A batch is screened against the box drawn round its own points. Building
  // that box one grid step short on any axis would drop a shell sitting on the
  // far face, and the points there would lose it entirely.
  //
  // One step has to be decisive for this to be visible at all. A shell only
  // dropped by a one-step shrink lies between cutoff-step and cutoff of the
  // box, where by construction it contributes about the shell tolerance, so a
  // fine grid hides the error inside the tolerance it is measured against.
  // The spacing is therefore set to 1.5x the largest cutoff radius, measured
  // from the basis rather than assumed, which puts a shell on a tile face
  // either fully in or fully out.
  //
  // Which grid point lands on a tile face depends on how the grid is
  // subdivided, which is not observable from here and moves with the batch
  // size, so the atom is swept over every grid point instead. Six points an
  // axis puts it on a low face, a far face and an interior point of every
  // axis under any tiling the batch size can produce.
  constexpr double shell_tol = 1e-10;

  double radius = 0.0;
  {
    Molecule probe;
    probe.emplace_back(AtomicNumber(8), 0.0, 0.0, 0.0);
    auto pb = make_ccpvdz(probe, SphericalType(true));
    for (auto& sh : pb) sh.set_shell_tolerance(shell_tol);
    for (const auto& sh : pb) radius = std::max(radius, sh.cutoff_radius());
  }
  REQUIRE(radius > 1.0);
  const double step = 1.5 * radius;

  CubeGrid grid;
  grid.origin = {0.0, 0.0, 0.0};
  grid.spacing = {step, step, step};
  grid.nx = 6;
  grid.ny = 6;
  grid.nz = 6;
  const int64_t npts = grid.num_points();
  const auto pts = grid.points();

  double worst_orb = 0.0, worst_rho = 0.0, max_rho = 0.0;
  for (int64_t site = 0; site < npts; ++site) {
    Molecule mol;
    mol.emplace_back(AtomicNumber(8), pts[3 * site + 0], pts[3 * site + 1],
                     pts[3 * site + 2]);

    auto basis = make_ccpvdz(mol, SphericalType(true));
    for (auto& sh : basis) sh.set_shell_tolerance(shell_tol);
    const int32_t nbf = basis.nbf();
    auto eval = make_evaluator(basis, shell_tol);
    const auto ao = reference_collocation(basis, npts, pts.data());

    const auto C = make_random_vector(static_cast<size_t>(nbf), 606u);
    std::vector<double> D(static_cast<size_t>(nbf) * nbf, 0.0);
    for (int32_t i = 0; i < nbf; ++i) D[static_cast<size_t>(i) * nbf + i] = 1.0;

    std::vector<double> orb(static_cast<size_t>(npts));
    std::vector<double> rho(static_cast<size_t>(npts));
    eval.eval_orbital(grid, C.data(), orb.data());
    eval.eval_density(grid, D.data(), nbf, rho.data());

    for (int64_t p = 0; p < npts; ++p) {
      const double* a = ao.data() + static_cast<size_t>(p) * nbf;
      double orb_ref = 0.0, rho_ref = 0.0;
      for (int32_t mu = 0; mu < nbf; ++mu) {
        orb_ref += C[static_cast<size_t>(mu)] * a[mu];
        rho_ref += a[mu] * a[mu];
      }
      max_rho = std::max(max_rho, rho_ref);
      worst_orb = std::max(worst_orb, std::fabs(orb[p] - orb_ref));
      worst_rho = std::max(worst_rho, std::fabs(rho[p] - rho_ref));
    }
  }

  CHECK(worst_orb < 1e-9);
  CHECK(worst_rho < 1e-9);
  // The atom has to register on the grid, or nothing is being proved about
  // whether its shells survived screening.
  CHECK(max_rho > 1e-2);
}

TEST_CASE("OrbitalEvaluator screens a tile spanning part of a long axis",
          "[orbital_evaluator]") {
  // A tile covers a sub-range of each axis, so on a grid far longer in z than
  // the tile is, the box has to close around that sub-range rather than the
  // whole axis. Three things have to line up for that to matter. The axis must
  // outlast a tile, so nz is far above any tile extent. The grid must span
  // well past a cutoff radius, taken here from the basis rather than
  // hard-coded. And the molecule must sit at the far end of z, so that a box
  // reaching further than its own points changes what survives screening.
  auto mol = make_water();

  constexpr double shell_tol = 1e-10;
  auto basis = make_ccpvdz(mol, SphericalType(true));
  for (auto& sh : basis) sh.set_shell_tolerance(shell_tol);
  const int32_t nbf = basis.nbf();

  std::vector<double> radii;
  radii.reserve(basis.size());
  for (const auto& sh : basis) radii.push_back(sh.cutoff_radius());
  const size_t mid = radii.size() / 2;
  std::nth_element(radii.begin(), radii.begin() + mid, radii.end());
  const double median_radius = radii[mid];

  const double z_extent = 4.0 * median_radius;

  CubeGrid grid;
  grid.nx = 1;
  grid.ny = 3;
  grid.nz = 3000;
  grid.spacing = {1.0, 1.0, z_extent / static_cast<double>(grid.nz)};
  // Offset in x so no point lands on a nucleus; the row runs from z_extent
  // below the molecule up to 1 Bohr short of it.
  grid.origin = {0.3, 0.0, -(z_extent + 1.0)};

  const int64_t npts = grid.num_points();
  const auto pts = grid.points();
  const auto ao_ref = reference_collocation(basis, npts, pts.data());

  const auto C = make_random_vector(static_cast<size_t>(nbf), 4099u);
  std::vector<double> D(static_cast<size_t>(nbf) * nbf, 0.0);
  for (int32_t i = 0; i < nbf; ++i) D[i * nbf + i] = 1.0;

  std::vector<double> orb(static_cast<size_t>(npts));
  std::vector<double> rho(static_cast<size_t>(npts));
  std::vector<double> orb_pts(static_cast<size_t>(npts));
  std::vector<double> rho_pts(static_cast<size_t>(npts));

  auto eval = make_evaluator(basis, shell_tol);
  eval.eval_orbital(grid, C.data(), orb.data());
  eval.eval_density(grid, D.data(), nbf, rho.data());
  // The pointer overload boxes the points it is handed instead of deriving a
  // box from tile indices, so it checks the index arithmetic independently.
  eval.eval_orbital(npts, pts.data(), C.data(), orb_pts.data());
  eval.eval_density(npts, pts.data(), D.data(), nbf, rho_pts.data());

  for (int64_t p = 0; p < npts; ++p) {
    double orb_ref = 0.0, rho_ref = 0.0;
    for (int32_t mu = 0; mu < nbf; ++mu) {
      const double a = ao_ref[static_cast<size_t>(p) * nbf + mu];
      orb_ref += C[static_cast<size_t>(mu)] * a;
      rho_ref += a * a;
    }
    CHECK(orb[p] == Approx(orb_ref).margin(1e-9));
    CHECK(rho[p] == Approx(rho_ref).margin(1e-9));
    CHECK(orb[p] == Approx(orb_pts[p]).margin(shell_tol));
    CHECK(rho[p] == Approx(rho_pts[p]).margin(shell_tol));
  }

  // Pins the span the tile box is being asked to resolve: the row ends on
  // the molecule and starts far enough away for shells to have died off.
  const size_t row_end = static_cast<size_t>(grid.nz) - 1;
  CHECK(rho[row_end] > 1e-2);
  CHECK(rho[0] < 1e-4 * rho[row_end]);
}

TEST_CASE("OrbitalEvaluator survives a thread-count change after construction",
          "[orbital_evaluator]") {
  // Regression guard: scratch must follow the thread count in force at
  // evaluation time, not at construction. Nothing about the batch
  // decomposition depends on the thread count, so the two runs screen
  // identically and have to agree bitwise.
  auto mol = make_water();
  auto basis = make_ccpvdz(mol, SphericalType(true));
  constexpr double shell_tol = 1e-12;
  for (auto& sh : basis) sh.set_shell_tolerance(shell_tol);
  const int32_t nbf = basis.nbf();

  auto grid = CubeGrid::from_molecule(mol, 16, 16, 16);
  const int64_t npts = grid.num_points();

  const auto C = make_random_vector(static_cast<size_t>(nbf), 5150u);
  std::vector<double> D(static_cast<size_t>(nbf) * nbf, 0.0);
  for (int32_t i = 0; i < nbf; ++i) D[i * nbf + i] = 1.0;

#ifdef _OPENMP
  const int saved_threads = omp_get_max_threads();
  // Construct while the thread count is 1, then raise it. Per-call scratch
  // must follow the thread count in force at evaluation time.
  omp_set_num_threads(1);
#endif
  auto eval = make_evaluator(basis, shell_tol);

  std::vector<double> orb_serial(static_cast<size_t>(npts));
  std::vector<double> rho_serial(static_cast<size_t>(npts));
  eval.eval_orbital(grid, C.data(), orb_serial.data());
  eval.eval_density(grid, D.data(), nbf, rho_serial.data());

#ifdef _OPENMP
  omp_set_num_threads(saved_threads > 1 ? saved_threads : 2);
#endif

  std::vector<double> orb_par(static_cast<size_t>(npts));
  std::vector<double> rho_par(static_cast<size_t>(npts));
  eval.eval_orbital(grid, C.data(), orb_par.data());
  eval.eval_density(grid, D.data(), nbf, rho_par.data());

#ifdef _OPENMP
  omp_set_num_threads(saved_threads);
#endif

  for (int64_t p = 0; p < npts; ++p) {
    CHECK(orb_par[p] == orb_serial[p]);
    CHECK(rho_par[p] == rho_serial[p]);
  }
}

TEST_CASE("OrbitalEvaluator is safe to invoke concurrently through const",
          "[orbital_evaluator]") {
  // The evaluator holds no mutable state, so a const invocation has to be safe
  // both from several threads at once and from inside an existing parallel
  // region. Neither mode crashed when scratch was shared between them; they
  // silently returned wrong numbers, which is why they are pinned here rather
  // than left to the thread-count test to catch.
  //
  // Nothing here changes the thread count, so every run derives the same batch
  // shape and screens identically. Agreement is therefore bitwise: any
  // difference at all means state leaked between calls.
  auto mol = make_water();
  constexpr double shell_tol = 1e-10;
  auto basis = make_ccpvdz(mol, SphericalType(true));
  for (auto& sh : basis) sh.set_shell_tolerance(shell_tol);
  const int32_t nbf = basis.nbf();
  auto eval = make_evaluator(basis, shell_tol);
  const OrbitalEvaluator& ceval = eval;

  auto grid = CubeGrid::from_molecule(mol, 12, 12, 12);
  const auto npts = static_cast<size_t>(grid.num_points());

  const auto C = make_random_vector(static_cast<size_t>(nbf), 271828u);
  std::vector<double> D(static_cast<size_t>(nbf) * nbf, 0.0);
  for (int32_t i = 0; i < nbf; ++i) D[static_cast<size_t>(i) * nbf + i] = 1.0;

  auto run = [&](double* orb, double* rho) {
    ceval.eval_orbital(grid, C.data(), orb);
    ceval.eval_density(grid, D.data(), nbf, rho);
  };

  std::vector<double> orb_ref(npts), rho_ref(npts);
  run(orb_ref.data(), rho_ref.data());
  double max_ref = 0.0;
  for (size_t p = 0; p < npts; ++p) max_ref = std::max(max_ref, rho_ref[p]);
  REQUIRE(max_ref > 1e-2);

  constexpr int nrunners = 3;
  std::vector<double> orb(npts * nrunners), rho(npts * nrunners);

  SECTION("concurrent std::threads on the same evaluator") {
    std::vector<std::thread> pool;
    for (int t = 0; t < nrunners; ++t)
      pool.emplace_back([&, t] {
        run(orb.data() + t * npts, rho.data() + t * npts);
      });
    for (auto& th : pool) th.join();
  }

#ifdef _OPENMP
  SECTION("evaluation from inside an enclosing parallel region") {
#pragma omp parallel for num_threads(nrunners) schedule(static, 1)
    for (int t = 0; t < nrunners; ++t)
      run(orb.data() + static_cast<size_t>(t) * npts,
          rho.data() + static_cast<size_t>(t) * npts);
  }
#endif

  size_t orb_diff = 0, rho_diff = 0;
  for (size_t t = 0; t < nrunners; ++t)
    for (size_t p = 0; p < npts; ++p) {
      if (orb[t * npts + p] != orb_ref[p]) ++orb_diff;
      if (rho[t * npts + p] != rho_ref[p]) ++rho_diff;
    }
  CHECK(orb_diff == 0);
  CHECK(rho_diff == 0);
}

TEST_CASE("OrbitalEvaluator results do not depend on the thread count",
          "[orbital_evaluator]") {
  // The batch decomposition is derived from the point count and the basis
  // alone, so every thread count screens against the same set of bounding
  // boxes and has to return the same bits. Deriving it from the thread count
  // instead is the easy mistake, and two things have to hold for that mistake
  // to be visible. Screening must be active, hence two centres 20 Bohr apart
  // with plenty of empty grid around them. And the batch size must be free to
  // move: below kMinBatches * kMinBatch points the minimum-batch-count bound
  // rounds away and the batch floor pins the batch whatever the thread count,
  // so the grid is deliberately larger than that.
  Molecule mol;
  mol.emplace_back(AtomicNumber(8), 0.0, 0.0, 0.0);
  mol.emplace_back(AtomicNumber(1), 20.0, 0.0, 0.0);
  auto basis = make_ccpvdz(mol, SphericalType(true));
  constexpr double shell_tol = 1e-10;
  for (auto& sh : basis) sh.set_shell_tolerance(shell_tol);
  const int32_t nbf = basis.nbf();
  auto eval = make_evaluator(basis, shell_tol);

  CubeGrid grid;
  grid.origin = {-4.0, -4.0, -4.0};
  grid.spacing = {0.2, 0.5, 0.5};
  grid.nx = 160;
  grid.ny = 32;
  grid.nz = 32;
  const auto npts = static_cast<size_t>(grid.num_points());
  REQUIRE(npts > 65536);

  const auto C = make_random_vector(static_cast<size_t>(nbf), 1618u);
  std::vector<double> D(static_cast<size_t>(nbf) * nbf, 0.0);
  for (int32_t i = 0; i < nbf; ++i) D[i * nbf + i] = 1.0;

  std::vector<double> orb_ref(npts), rho_ref(npts);
  std::vector<double> orb(npts), rho(npts);

#ifdef _OPENMP
  const int saved_threads = omp_get_max_threads();
  omp_set_num_threads(1);
#endif
  eval.eval_orbital(grid, C.data(), orb_ref.data());
  eval.eval_density(grid, D.data(), nbf, rho_ref.data());

  // Screening has to be doing something, or the invariance is vacuous.
  double max_rho = 0.0;
  for (size_t p = 0; p < npts; ++p) max_rho = std::max(max_rho, rho_ref[p]);
  REQUIRE(max_rho > 1e-2);

  size_t orb_diff = 0, rho_diff = 0;
  for (int nthreads : {2, 3, 5, 8, 13, 16, 32}) {
#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#else
    (void)nthreads;
#endif
    eval.eval_orbital(grid, C.data(), orb.data());
    eval.eval_density(grid, D.data(), nbf, rho.data());
    for (size_t p = 0; p < npts; ++p) {
      if (orb[p] != orb_ref[p]) ++orb_diff;
      if (rho[p] != rho_ref[p]) ++rho_diff;
    }
  }
#ifdef _OPENMP
  omp_set_num_threads(saved_threads);
#endif

  CHECK(orb_diff == 0);
  CHECK(rho_diff == 0);
}

TEST_CASE("OrbitalEvaluator screening error scales with the shell tolerance",
          "[orbital_evaluator]") {
  // Pins the guidance on the class: orbital error is bounded by the shell
  // tolerance, density error by its square, since dropping a shell with
  // |phi| < t perturbs sum_uv D phi phi by ~t^2.
  Molecule mol;
  mol.emplace_back(AtomicNumber(8), 0.0, 0.0, 0.0);
  mol.emplace_back(AtomicNumber(1), 20.0, 0.0, 0.0);

  constexpr double loose_tol = 1e-6;
  auto basis = make_ccpvdz(mol, SphericalType(true));
  const int32_t nbf = basis.nbf();

  std::vector<double> radii_before;
  for (const auto& sh : basis) radii_before.push_back(sh.cutoff_radius());

  auto eval_tight = make_evaluator(basis, 1e-14);
  auto eval_loose = make_evaluator(basis, loose_tol);

  // Retuning happens on each evaluator's own copy; a basis shared with an SCF
  // setup must come back unchanged.
  for (size_t i = 0; i < basis.size(); ++i)
    CHECK(basis[i].cutoff_radius() == radii_before[i]);

  CubeGrid grid;
  grid.origin = {-4.0, -4.0, -4.0};
  grid.spacing = {0.8, 2.0, 2.0};
  grid.nx = 40;
  grid.ny = 4;
  grid.nz = 4;
  const int64_t npts = grid.num_points();

  const auto C = make_random_vector(static_cast<size_t>(nbf), 31337u);
  std::vector<double> D(static_cast<size_t>(nbf) * nbf, 0.0);
  for (int32_t i = 0; i < nbf; ++i) D[i * nbf + i] = 1.0;

  std::vector<double> orb_tight(npts), orb_loose(npts);
  std::vector<double> rho_tight(npts), rho_loose(npts);
  eval_tight.eval_orbital(grid, C.data(), orb_tight.data());
  eval_loose.eval_orbital(grid, C.data(), orb_loose.data());
  eval_tight.eval_density(grid, D.data(), nbf, rho_tight.data());
  eval_loose.eval_density(grid, D.data(), nbf, rho_loose.data());

  double orb_err = 0.0, rho_err = 0.0;
  for (int64_t p = 0; p < npts; ++p) {
    orb_err = std::max(orb_err, std::fabs(orb_loose[p] - orb_tight[p]));
    rho_err = std::max(rho_err, std::fabs(rho_loose[p] - rho_tight[p]));
  }

  INFO("orbital error " << orb_err << ", density error " << rho_err);
  CHECK(orb_err > 0.0);  // screening is genuinely active at the loose tolerance
  CHECK(orb_err <= 10.0 * loose_tol);
  CHECK(rho_err <= 100.0 * loose_tol * loose_tol);
}

TEST_CASE("CubeGrid construction and grid-points layout", "[cube]") {
  auto mol = make_water();
  CubeGrid g = CubeGrid::from_molecule(mol, /*nx=*/16, /*ny=*/12, /*nz=*/8,
                                       /*margin=*/2.0);
  REQUIRE(g.num_points() == 16 * 12 * 8);

  auto pts = g.points();
  REQUIRE(pts.size() == static_cast<size_t>(g.num_points()) * 3);

  // Check first and last grid points.
  CHECK(pts[0] == Approx(g.origin[0]));
  CHECK(pts[1] == Approx(g.origin[1]));
  CHECK(pts[2] == Approx(g.origin[2]));

  const size_t last = static_cast<size_t>(g.num_points()) - 1;
  CHECK(pts[3 * last + 0] ==
        Approx(g.origin[0] + g.spacing[0] * (g.nx - 1)));
  CHECK(pts[3 * last + 1] ==
        Approx(g.origin[1] + g.spacing[1] * (g.ny - 1)));
  CHECK(pts[3 * last + 2] ==
        Approx(g.origin[2] + g.spacing[2] * (g.nz - 1)));

  // Check ordering: iz varies fastest. Point (1, 0, 0) should be at offset
  // ny*nz, point (0, 1, 0) at nz, point (0, 0, 1) at 1.
  const int64_t off_x = g.ny * g.nz;
  const int64_t off_y = g.nz;
  CHECK(pts[3 * static_cast<size_t>(off_x) + 0] ==
        Approx(g.origin[0] + g.spacing[0]));
  CHECK(pts[3 * static_cast<size_t>(off_y) + 1] ==
        Approx(g.origin[1] + g.spacing[1]));
  CHECK(pts[3 * 1 + 2] == Approx(g.origin[2] + g.spacing[2]));
}

TEST_CASE("CubeGrid rejects degenerate specifications", "[cube]") {
  auto mol = make_water();
  CHECK_THROWS(CubeGrid::from_molecule(Molecule{}, 4, 4, 4));
  CHECK_THROWS(CubeGrid::from_molecule(mol, 0, 4, 4));
  CHECK_THROWS(CubeGrid::from_molecule(mol, 4, 0, 4));
  CHECK_THROWS(CubeGrid::from_molecule(mol, 4, 4, 0));

  const auto grid = CubeGrid::from_molecule(mol, 3, 4, 5);
  const auto pts = grid.points();
  REQUIRE(pts.size() == static_cast<size_t>(grid.num_points()) * 3);

  std::vector<double> buf(pts.size(), -1.0);
  grid.points_into(buf.data());
  for (size_t i = 0; i < pts.size(); ++i) CHECK(buf[i] == pts[i]);
}

TEST_CASE("write_cube rejects bad input and defaults its comment", "[cube]") {
#ifdef GAUXC_HAS_MPI
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_rank) return;  // File I/O; only run on root rank
#endif
  auto mol = make_water();
  const auto grid = CubeGrid::from_molecule(mol, 2, 2, 2);
  std::vector<double> field(static_cast<size_t>(grid.num_points()), 1.0);
  const auto path = make_temp_path(".cube");

  CHECK_THROWS(write_cube(path, mol, grid, nullptr));

  CubeGrid empty;
  empty.nx = 0;
  CHECK_THROWS(write_cube(path, mol, empty, field.data()));
  CHECK_THROWS(
      write_cube("/nonexistent-gauxc-dir/out.cube", mol, grid, field.data()));

#ifdef GAUXC_HAS_HDF5
  const auto h5 = make_temp_path(".h5");
  CHECK_THROWS(write_cube_hdf5(h5, mol, grid, nullptr));
  CHECK_THROWS(write_cube_hdf5(h5, mol, empty, field.data()));
#endif

  // An omitted comment falls back to a fixed first line.
  write_cube(path, mol, grid, field.data());
  std::ifstream in(path);
  REQUIRE(in.good());
  std::string line;
  std::getline(in, line);
  CHECK(line == "GauXC cube file");
  std::getline(in, line);
  CHECK(line == "Generated by GauXC");
}

TEST_CASE("write_cube round-trips header and field data", "[cube]") {
#ifdef GAUXC_HAS_MPI
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_rank) return;  // File I/O; only run on root rank
#endif
  auto mol = make_water();
  CubeGrid grid = CubeGrid::from_molecule(mol, /*nx=*/5, /*ny=*/4, /*nz=*/7,
                                          /*margin=*/2.0);
  std::vector<double> field(static_cast<size_t>(grid.num_points()));
  for (int64_t i = 0; i < grid.num_points(); ++i) {
    // Mix of magnitudes to exercise the formatter (negative, near-zero, etc.)
    field[static_cast<size_t>(i)] = std::sin(0.13 * static_cast<double>(i)) *
                                    std::pow(10.0, (i % 5) - 2);
  }

  const std::string path = make_temp_path(".cube");
  write_cube(path, mol, grid, field.data(), "Test cube");

  // Parse back.
  std::ifstream in(path);
  REQUIRE(in.is_open());
  std::string l1, l2;
  std::getline(in, l1);
  std::getline(in, l2);
  CHECK(l1 == "Test cube");
  CHECK(l2 == "Generated by GauXC");

  long long natoms_read;
  double ox, oy, oz;
  in >> natoms_read >> ox >> oy >> oz;
  CHECK(natoms_read == static_cast<long long>(mol.size()));
  CHECK(ox == Approx(grid.origin[0]));
  CHECK(oy == Approx(grid.origin[1]));
  CHECK(oz == Approx(grid.origin[2]));

  // 3 axis lines.
  for (int axis = 0; axis < 3; ++axis) {
    long long n;
    double a, b, c;
    in >> n >> a >> b >> c;
    if (axis == 0) {
      CHECK(n == grid.nx);
      CHECK(a == Approx(grid.spacing[0]));
      CHECK(b == Approx(0.0));
      CHECK(c == Approx(0.0));
    } else if (axis == 1) {
      CHECK(n == grid.ny);
      CHECK(a == Approx(0.0));
      CHECK(b == Approx(grid.spacing[1]));
      CHECK(c == Approx(0.0));
    } else {
      CHECK(n == grid.nz);
      CHECK(a == Approx(0.0));
      CHECK(b == Approx(0.0));
      CHECK(c == Approx(grid.spacing[2]));
    }
  }

  // Atoms.
  for (size_t i = 0; i < mol.size(); ++i) {
    long long Z;
    double q, x, y, z;
    in >> Z >> q >> x >> y >> z;
    CHECK(Z == mol[i].Z.get());
    CHECK(q == Approx(0.0));
    CHECK(x == Approx(mol[i].x));
    CHECK(y == Approx(mol[i].y));
    CHECK(z == Approx(mol[i].z));
  }

  // Field. Read all remaining whitespace-separated doubles and compare.
  std::vector<double> field_read;
  field_read.reserve(static_cast<size_t>(grid.num_points()));
  double v;
  while (in >> v) field_read.push_back(v);

  REQUIRE(field_read.size() == field.size());
  // %13.5E gives 5 significant digits → relative tolerance ~1e-5 for the
  // round-trip.
  for (size_t i = 0; i < field.size(); ++i) {
    CHECK(field_read[i] == Approx(field[i]).epsilon(1e-4).margin(1e-30));
  }
  in.close();

  // Line structure: each (ix,iy) row spans ceil(nz/6) lines, full lines carry
  // six 13-char fields and the row's last line carries the remainder. This
  // pins the exact per-row byte count that write_cube relies on to place rows
  // without a compaction pass.
  {
    std::ifstream lin(path);
    REQUIRE(lin.is_open());
    std::string line;
    for (size_t i = 0; i < 6 + mol.size(); ++i) std::getline(lin, line);

    const int64_t lines_per_row = (grid.nz + 5) / 6;
    for (int64_t row = 0; row < grid.nx * grid.ny; ++row) {
      for (int64_t l = 0; l < lines_per_row; ++l) {
        REQUIRE(static_cast<bool>(std::getline(lin, line)));
        const int64_t nvals = (l + 1 == lines_per_row) ? (grid.nz - l * 6) : 6;
        CHECK(line.size() == static_cast<size_t>(13 * nvals));
      }
    }
    CHECK_FALSE(static_cast<bool>(std::getline(lin, line)));
  }

  std::remove(path.c_str());
}

namespace {

/// Write `field` as a single-row cube file and return the raw data block.
std::string cube_data_block(const Molecule& mol,
                            const std::vector<double>& field) {
  CubeGrid grid;
  grid.origin = {0.0, 0.0, 0.0};
  grid.spacing = {0.1, 0.1, 0.1};
  grid.nx = 1;
  grid.ny = 1;
  grid.nz = static_cast<int64_t>(field.size());

  const std::string path = make_temp_path(".cube");
  write_cube(path, mol, grid, field.data(), "fmt");

  std::ifstream in(path);
  REQUIRE(in.is_open());
  // Skip header (2 comments + 1 natoms line + 3 axis lines + natoms atoms).
  std::string skip;
  for (int i = 0; i < 6; ++i) std::getline(in, skip);
  for (size_t i = 0; i < mol.size(); ++i) std::getline(in, skip);

  std::string block((std::istreambuf_iterator<char>(in)),
                    std::istreambuf_iterator<char>());
  in.close();
  std::remove(path.c_str());
  return block;
}

}  // namespace

TEST_CASE("write_cube agrees with snprintf %13.5E formatting", "[cube]") {
#ifdef GAUXC_HAS_MPI
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_rank) return;  // File I/O; only run on root rank
#endif
  // Spot-check the custom formatter against snprintf for a battery of values
  // by writing a tiny cube file and parsing it back. This is a stronger check
  // than the round-trip above since we compare the byte stream.
  auto mol = make_water();

  const double qnan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  // The last five exercise the snprintf deferral: subnormals (where scaling
  // the mantissa would overflow), a 3-digit exponent either side of zero, and
  // a mantissa that carries from E+99 up into a 3-digit exponent.
  std::vector<double> field = {0.0,
                               -0.0,
                               1.23456e-10,
                               -9.99995e-1,
                               1.0e+99,
                               -1.0e+99,
                               3.14159265358979,
                               -2.71828,
                               1.0,
                               -1.0,
                               1e-300,
                               1.234e+05,
                               qnan,
                               -qnan,
                               inf,
                               -inf,
                               std::numeric_limits<double>::denorm_min(),
                               -1e-310,
                               1.0e+300,
                               -1.0e-305,
                               9.9999999e+99};

  const std::string data_block = cube_data_block(mol, field);

  std::ostringstream expected;
  for (size_t i = 0; i < field.size(); ++i) {
    char buf[32];
    // C99 leaves the sign of a printed NaN implementation-defined: glibc emits
    // it, BSD libc does not, so snprintf is not a portable reference for that
    // one value. The formatter deliberately emits it, which is pinned here.
    if (std::isnan(field[i]))
      std::snprintf(buf, sizeof(buf), "%13s",
                    std::signbit(field[i]) ? "-NAN" : "NAN");
    else
      std::snprintf(buf, sizeof(buf), "%13.5E", field[i]);
    expected << buf;
    if ((i + 1) % 6 == 0 || (i + 1) == field.size()) expected << '\n';
  }

  CHECK(data_block == expected.str());
}

TEST_CASE("write_cube spans multiple output chunks", "[cube]") {
#ifdef GAUXC_HAS_MPI
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_rank) return;  // File I/O; only run on root rank
#endif
  // Sized to exceed write_cube's internal staging-buffer target so the chunk
  // loop runs more than once, with a partial final chunk.
  auto mol = make_water();
  CubeGrid grid;
  grid.origin = {0.0, 0.0, 0.0};
  grid.spacing = {0.1, 0.1, 0.1};
  grid.nx = 2;
  grid.ny = 100;
  grid.nz = 4096;
  const int64_t npts = grid.num_points();

  std::vector<double> field(static_cast<size_t>(npts));
  for (int64_t i = 0; i < npts; ++i)
    field[static_cast<size_t>(i)] = std::sin(1e-3 * static_cast<double>(i));

  const std::string path = make_temp_path(".cube");
  write_cube(path, mol, grid, field.data(), "chunked");

  std::ifstream in(path, std::ios::binary);
  REQUIRE(in.is_open());
  std::string skip;
  for (size_t i = 0; i < 6 + mol.size(); ++i) std::getline(in, skip);

  std::string data_block((std::istreambuf_iterator<char>(in)),
                         std::istreambuf_iterator<char>());
  in.close();

  const int64_t bytes_per_row = grid.nz * 13 + (grid.nz + 5) / 6;
  const int64_t n_rows = grid.nx * grid.ny;
  REQUIRE(data_block.size() == static_cast<size_t>(bytes_per_row * n_rows));

  // Full byte comparison catches any misplacement at a chunk seam.
  std::string expected;
  expected.reserve(data_block.size());
  char buf[32];
  for (int64_t row = 0; row < n_rows; ++row) {
    for (int64_t iz = 0; iz < grid.nz; ++iz) {
      std::snprintf(buf, sizeof(buf), "%13.5E",
                    field[static_cast<size_t>(row * grid.nz + iz)]);
      expected += buf;
      if ((iz + 1) % 6 == 0 || iz + 1 == grid.nz) expected += '\n';
    }
  }
  CHECK(data_block == expected);

  std::remove(path.c_str());
}

TEST_CASE("write_cube formatter stays within one last digit at a rounding "
          "boundary",
          "[cube]") {
#ifdef GAUXC_HAS_MPI
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_rank) return;  // File I/O; only run on root rank
#endif
  // For values sitting within a few ulp of a rounding boundary in the 6th
  // significant digit, scaling the mantissa can tip it across the boundary,
  // so the formatter and glibc may pick different last digits. Both stay
  // within one unit of it; this pins that bound rather than byte equality.
  auto mol = make_water();
  const std::vector<double> field = {123456.5, 1.234575, -123456.5, -1.234575,
                                     9.9999949999999998642e-98, 0.0};

  const std::string data_block = cube_data_block(mol, field);
  REQUIRE(data_block.size() >= field.size() * 13);

  for (size_t i = 0; i < field.size(); ++i) {
    const std::string tok = data_block.substr(i * 13, 13);
    const double got = std::stod(tok);
    const double v = field[i];
    const double last_digit =
        v == 0.0 ? 1.0
                 : std::pow(10.0, std::floor(std::log10(std::fabs(v))) - 5.0);
    INFO("value " << v << " formatted as '" << tok << "'");
    CHECK(std::fabs(got - v) <= 0.5000001 * last_digit);
  }
}

#ifdef GAUXC_HAS_HDF5
TEST_CASE("write_cube_hdf5 round-trip", "[cube]") {
#ifdef GAUXC_HAS_MPI
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_rank) return;  // File I/O; only run on root rank
#endif
  auto mol = make_water();
  auto grid = CubeGrid::from_molecule(mol, 4, 5, 6);

  // Fill a small field with known values.
  const int64_t npts = grid.num_points();
  std::vector<double> field(npts);
  for (int64_t i = 0; i < npts; ++i) field[i] = 0.01 * i - 0.5;

  const std::string path = make_temp_path(".h5");
  write_cube_hdf5(path, mol, grid, field.data(), "test cube");

  // Read back and verify.
  HighFive::File file(path, HighFive::File::ReadOnly);

  // Field shape and values.
  auto ds = file.getDataSet("field");
  auto dims = ds.getDimensions();
  REQUIRE(dims.size() == 3);
  CHECK(dims[0] == static_cast<size_t>(grid.nx));
  CHECK(dims[1] == static_cast<size_t>(grid.ny));
  CHECK(dims[2] == static_cast<size_t>(grid.nz));

  std::vector<double> read_field(npts);
  ds.read(read_field.data());
  for (int64_t i = 0; i < npts; ++i) {
    CHECK(read_field[i] == Approx(field[i]).epsilon(1e-14));
  }

  // Comment attribute.
  std::string cmt;
  ds.getAttribute("comment").read(cmt);
  CHECK(cmt == "test cube");

  // Grid metadata.
  auto grp_grid = file.getGroup("grid");
  std::vector<double> origin, spacing;
  std::vector<int64_t> shape;
  grp_grid.getDataSet("origin").read(origin);
  grp_grid.getDataSet("spacing").read(spacing);
  grp_grid.getDataSet("shape").read(shape);
  REQUIRE(origin.size() == 3);
  REQUIRE(spacing.size() == 3);
  REQUIRE(shape.size() == 3);
  for (int k = 0; k < 3; ++k) {
    CHECK(origin[k] == Approx(grid.origin[k]).epsilon(1e-14));
    CHECK(spacing[k] == Approx(grid.spacing[k]).epsilon(1e-14));
  }
  CHECK(shape[0] == grid.nx);
  CHECK(shape[1] == grid.ny);
  CHECK(shape[2] == grid.nz);

  // Atoms.
  auto grp_atoms = file.getGroup("atoms");
  std::vector<int64_t> Z;
  grp_atoms.getDataSet("Z").read(Z);
  REQUIRE(Z.size() == mol.size());
  for (size_t i = 0; i < mol.size(); ++i) {
    CHECK(Z[i] == static_cast<int64_t>(mol[i].Z.get()));
  }

  std::vector<double> coords(mol.size() * 3);
  grp_atoms.getDataSet("coords").read(coords.data());
  for (size_t i = 0; i < mol.size(); ++i) {
    CHECK(coords[3 * i + 0] == Approx(mol[i].x).epsilon(1e-14));
    CHECK(coords[3 * i + 1] == Approx(mol[i].y).epsilon(1e-14));
    CHECK(coords[3 * i + 2] == Approx(mol[i].z).epsilon(1e-14));
  }

  std::remove(path.c_str());
}
#endif
