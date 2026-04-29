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
#include "ut_common.hpp"
#include "catch2/catch.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <random>
#include <sstream>
#include <vector>

#include <unistd.h>

#include <gauxc/external/cube.hpp>
#include <gauxc/orbital_evaluator.hpp>
#include <gauxc/xc_integrator/local_work_driver.hpp>

#include "standards.hpp"

// Reference implementation lives in src/, behind the public header surface.
#include "xc_integrator/local_work_driver/host/local_host_work_driver.hpp"

using namespace GauXC;

namespace {

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
  std::vector<double> ao_ref(static_cast<size_t>(nbf) * npts);
  {
    auto drv = LocalWorkDriverFactory::make_local_work_driver(
        ExecutionSpace::Host, "Reference");
    auto* host_drv = dynamic_cast<LocalHostWorkDriver*>(drv.get());
    REQUIRE(host_drv != nullptr);
    std::vector<int32_t> shell_list(basis.size());
    for (size_t i = 0; i < shell_list.size(); ++i)
      shell_list[i] = static_cast<int32_t>(i);
    host_drv->eval_collocation(static_cast<size_t>(npts),
                               static_cast<size_t>(basis.nshells()),
                               static_cast<size_t>(nbf), pts.data(), basis,
                               shell_list.data(), ao_ref.data());
  }

  OrbitalEvaluator eval(basis);
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
    std::vector<double> C(static_cast<size_t>(nbf) * nmo);
    std::mt19937 gen(99);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    for (auto& v : C) v = u(gen);

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
    std::vector<double> c(nbf);
    std::mt19937 gen(7);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    for (auto& v : c) v = u(gen);

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

TEST_CASE("write_cube round-trips header and field data", "[cube]") {
  auto mol = make_water();
  CubeGrid grid = CubeGrid::from_molecule(mol, /*nx=*/5, /*ny=*/4, /*nz=*/7,
                                          /*margin=*/2.0);
  std::vector<double> field(static_cast<size_t>(grid.num_points()));
  for (int64_t i = 0; i < grid.num_points(); ++i) {
    // Mix of magnitudes to exercise the formatter (negative, near-zero, etc.)
    field[static_cast<size_t>(i)] = std::sin(0.13 * static_cast<double>(i)) *
                                    std::pow(10.0, (i % 5) - 2);
  }

  // Use a tmp path under /tmp.
  const std::string path = std::string("/tmp/gauxc_cube_test_") +
                           std::to_string(::getpid()) + ".cube";
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

  std::remove(path.c_str());
}

TEST_CASE("write_cube agrees with snprintf %13.5E formatting", "[cube]") {
  // Spot-check the custom formatter against snprintf for a battery of values
  // by writing a tiny cube file and parsing it back. This is a stronger check
  // than the round-trip above since we compare the byte stream.
  auto mol = make_water();
  CubeGrid grid;
  grid.origin = {0.0, 0.0, 0.0};
  grid.spacing = {0.1, 0.1, 0.1};
  grid.nx = 1;
  grid.ny = 1;
  grid.nz = 12;

  std::vector<double> field = {0.0,        -0.0,    1.23456e-10, -9.99995e-1,
                               1.0e+99,    -1.0e+99,  3.14159265358979,
                               -2.71828,   1.0,      -1.0,       1e-300,
                               1.234e+05};

  const std::string path = std::string("/tmp/gauxc_cube_format_") +
                           std::to_string(::getpid()) + ".cube";
  write_cube(path, mol, grid, field.data(), "fmt");

  std::ifstream in(path);
  REQUIRE(in.is_open());
  // Skip header (2 comments + 4 grid lines + natoms atom lines).
  std::string skip;
  for (int i = 0; i < 2; ++i) std::getline(in, skip);
  std::getline(in, skip);  // natoms line
  for (int i = 0; i < 3; ++i) std::getline(in, skip);  // 3 axis lines
  for (size_t i = 0; i < mol.size(); ++i) std::getline(in, skip);

  // Read the field-block as raw text and compare against snprintf
  // line-by-line. Six values per line, then row terminator after 12 (single
  // (ix, iy) row here so no early newline).
  std::string data_block((std::istreambuf_iterator<char>(in)),
                         std::istreambuf_iterator<char>());

  std::ostringstream expected;
  for (size_t i = 0; i < field.size(); ++i) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%13.5E", field[i]);
    expected << buf;
    if ((i + 1) % 6 == 0 || (i + 1) == field.size()) expected << '\n';
  }

  CHECK(data_block == expected.str());

  std::remove(path.c_str());
}
