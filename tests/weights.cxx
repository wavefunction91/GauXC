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
#include <gauxc/molgrid.hpp>
#include <gauxc/molgrid/defaults.hpp>
#include <gauxc/basisset.hpp>
#include <gauxc/load_balancer.hpp>
#include <gauxc/util/div_ceil.hpp>
#include <cmath>
#include <fstream>
#include <string>

#include "weights_generate.hpp"
#include "weights_host.hpp"
#include "weights_cuda.hpp"
#include "weights_hip.hpp"

//#define GENERATE_TESTS
TEST_CASE( "Partition Weights", "[weights]" ) {

  auto rt = RuntimeEnvironment(GAUXC_MPI_CODE(MPI_COMM_WORLD));
#ifdef GENERATE_TESTS
  if(rt.comm_size() > 1) return;
#endif

  Molecule mol = make_benzene();

#ifdef GENERATE_TESTS
  BasisSet<double> basis = make_631Gd( mol, SphericalType(true) );
  for( auto& sh : basis ) sh.set_shell_tolerance( 1e-6 );

  generate_weights_data( mol, basis, "benzene_weights_becke.hdf5", XCWeightAlg::Becke );  
  generate_weights_data( mol, basis, "benzene_weights_ssf.hdf5", XCWeightAlg::SSF );  
  generate_weights_data( mol, basis, "benzene_weights_lko.hdf5", XCWeightAlg::LKO );
  generate_weights_data( mol, basis, "benzene_weights_hirshfeld.hdf5", XCWeightAlg::Hirshfeld );
  return;
#else


#ifdef GAUXC_HAS_HOST
  SECTION("Becke") {
  std::string ref_file = GAUXC_REF_DATA_PATH "/benzene_weights_becke.hdf5";
  test_host_weights( ref_file, XCWeightAlg::Becke );
  }
  SECTION("LKO") {
  std::string ref_file = GAUXC_REF_DATA_PATH "/benzene_weights_lko.hdf5";
  test_host_weights( ref_file, XCWeightAlg::LKO );
  }
  SECTION("Hirshfeld") {
  std::string ref_file = GAUXC_REF_DATA_PATH "/benzene_weights_hirshfeld.hdf5";
  test_host_weights( ref_file, XCWeightAlg::Hirshfeld );
  }
#endif


  SECTION("SSF") {

  std::string ref_file = GAUXC_REF_DATA_PATH "/benzene_weights_ssf.hdf5";

#ifdef GAUXC_HAS_HOST
  SECTION( "Host Weights" ) {
    test_host_weights( ref_file, XCWeightAlg::SSF );
  }
#endif

#ifdef GAUXC_HAS_DEVICE
  SECTION( "Device Weights" ) {
#ifdef GAUXC_HAS_CUDA
    test_cuda_weights( ref_file );
#elif defined(GAUXC_HAS_HIP)
    test_hip_weights( ref_file );
#endif
  }
#endif
#endif

  }
}

#ifdef GAUXC_HAS_HOST
namespace {

double hirshfeld_test_density(
  const Atom& atom,
  const std::array<double,3>& point
) {

  constexpr double pi = 3.141592653589793238462643383279502884;
  const double radius = default_atomic_radius(atom.Z);
  const double alpha  = 0.5 / (radius * radius);

  const double da_x = point[0] - atom.x;
  const double da_y = point[1] - atom.y;
  const double da_z = point[2] - atom.z;
  const double r2 = da_x*da_x + da_y*da_y + da_z*da_z;

  return atom.Z.get() * std::pow(alpha / pi, 1.5) * std::exp(-alpha * r2);
}

double hirshfeld_test_factor(
  const Molecule& mol,
  size_t parent,
  const std::array<double,3>& point
) {

  double sum = 0.;
  for( const auto& atom : mol ) sum += hirshfeld_test_density(atom, point);

  return (sum > 0.) ? hirshfeld_test_density(mol[parent], point) / sum : 0.;
}

Molecule translated_molecule(
  const Molecule& mol,
  const std::array<double,3>& shift
) {

  Molecule translated;
  translated.reserve(mol.size());
  for( const auto& atom : mol ) {
    translated.emplace_back(
      atom.Z, atom.x + shift[0], atom.y + shift[1], atom.z + shift[2]
    );
  }

  return translated;
}

std::vector<std::array<double,3>> translated_points(
  const std::vector<std::array<double,3>>& points,
  const std::array<double,3>& shift
) {

  std::vector<std::array<double,3>> translated = points;
  for( auto& point : translated ) {
    point[0] += shift[0];
    point[1] += shift[1];
    point[2] += shift[2];
  }

  return translated;
}

}

TEST_CASE( "Hirshfeld Partition Weights", "[weights]" ) {

  SECTION("H2 symmetric Gaussian pro-atom partition") {
    Molecule mol;
    mol.emplace_back( AtomicNumber(1), 0.0, 0.0, 0.0 );
    mol.emplace_back( AtomicNumber(1), 1.4, 0.0, 0.0 );
    MolMeta meta(mol);

    std::vector<XCTask> tasks(1);
    tasks[0].iParent = 0;
    tasks[0].points = {{
      {0.0, 0.0, 0.0},
      {0.7, 0.0, 0.0},
      {1.4, 0.0, 0.0}
    }};
    tasks[0].weights = {1.0, 1.0, 1.0};

    reference_hirshfeld_weights_host( mol, meta, tasks.begin(), tasks.end() );

    CHECK( tasks[0].weights[0] == Approx(0.987761411136141) );
    CHECK( tasks[0].weights[1] == Approx(0.5) );
    CHECK( tasks[0].weights[2] == Approx(0.0122385888638589) );
  }

  SECTION("CO heteronuclear Gaussian pro-atom partition") {
    Molecule mol;
    mol.emplace_back( AtomicNumber(6), 0.0, 0.0, 0.0 );
    mol.emplace_back( AtomicNumber(8), 2.2, 0.0, 0.0 );
    MolMeta meta(mol);

    std::vector<XCTask> tasks(2);
    tasks[0].iParent = 0;
    tasks[0].points = {{ {1.1, 0.0, 0.0} }};
    tasks[0].weights = {1.0};
    tasks[1].iParent = 1;
    tasks[1].points = tasks[0].points;
    tasks[1].weights = {1.0};

    reference_hirshfeld_weights_host( mol, meta, tasks.begin(), tasks.end() );

    CHECK( tasks[0].weights[0] == Approx(0.348581523757377) );
    CHECK( tasks[1].weights[0] == Approx(0.651418476242623) );
    CHECK( tasks[0].weights[0] + tasks[1].weights[0] == Approx(1.0) );
  }

  SECTION("single atom preserves quadrature weights") {
    Molecule mol;
    mol.emplace_back( AtomicNumber(8), -0.3, 0.4, 1.1 );
    MolMeta meta(mol);

    std::vector<XCTask> tasks(1);
    tasks[0].iParent = 0;
    tasks[0].points = {{
      {-0.3, 0.4, 1.1},
      { 0.2, 0.9, 1.7},
      {-1.5, 0.2, 0.1}
    }};
    tasks[0].weights = {0.25, 2.0, 7.5};
    const auto weights_ref = tasks[0].weights;

    reference_hirshfeld_weights_host( mol, meta, tasks.begin(), tasks.end() );

    REQUIRE( tasks[0].weights.size() == weights_ref.size() );
    for( size_t i = 0; i < weights_ref.size(); ++i ) {
      CHECK( tasks[0].weights[i] == Approx(weights_ref[i]) );
    }
  }

  SECTION("partition factors conserve each original point weight") {
    Molecule mol;
    mol.emplace_back( AtomicNumber(6), -0.7,  0.0,  0.2 );
    mol.emplace_back( AtomicNumber(1),  1.2, -0.4,  0.5 );
    mol.emplace_back( AtomicNumber(8),  0.3,  1.1, -0.3 );
    mol.emplace_back( AtomicNumber(7), -1.4, -0.8,  0.9 );
    MolMeta meta(mol);

    const std::vector<std::array<double,3>> points = {{
      {-0.1,  0.2,  0.0},
      { 0.9, -0.2,  0.4},
      {-1.0, -0.6,  1.3},
      { 1.8,  1.5, -0.9}
    }};
    const std::vector<double> point_weights = {0.25, 1.0, 3.5, 10.0};

    std::vector<XCTask> tasks(mol.natoms());
    for( size_t parent = 0; parent < mol.natoms(); ++parent ) {
      tasks[parent].iParent = static_cast<int32_t>(parent);
      tasks[parent].points = points;
      tasks[parent].weights = point_weights;
    }

    reference_hirshfeld_weights_host( mol, meta, tasks.begin(), tasks.end() );

    for( size_t ipt = 0; ipt < points.size(); ++ipt ) {
      double partition_sum = 0.;
      for( size_t parent = 0; parent < mol.natoms(); ++parent ) {
        const double factor = hirshfeld_test_factor(mol, parent, points[ipt]);
        CHECK( tasks[parent].weights[ipt] == Approx(point_weights[ipt] * factor) );
        partition_sum += tasks[parent].weights[ipt];
      }
      CHECK( partition_sum == Approx(point_weights[ipt]) );
    }
  }

  SECTION("rigid translation leaves partition factors unchanged") {
    Molecule mol;
    mol.emplace_back( AtomicNumber(6), -0.4, 0.0, 0.1 );
    mol.emplace_back( AtomicNumber(8),  1.6, 0.3, 0.0 );
    mol.emplace_back( AtomicNumber(1),  0.2, 1.5, 0.8 );

    const std::vector<std::array<double,3>> points = {{
      { 0.1, 0.1, 0.2},
      { 0.8, 0.4, 0.1},
      {-0.9, 1.2, 0.5}
    }};
    const std::vector<double> point_weights = {1.0, 2.0, 4.0};
    const std::array<double,3> shift = {8.0, -3.5, 2.25};

    auto shifted_mol = translated_molecule(mol, shift);
    auto shifted_points = translated_points(points, shift);

    MolMeta meta(mol);
    MolMeta shifted_meta(shifted_mol);

    std::vector<XCTask> tasks(1), shifted_tasks(1);
    tasks[0].iParent = 1;
    tasks[0].points = points;
    tasks[0].weights = point_weights;
    shifted_tasks[0].iParent = 1;
    shifted_tasks[0].points = shifted_points;
    shifted_tasks[0].weights = point_weights;

    reference_hirshfeld_weights_host( mol, meta, tasks.begin(), tasks.end() );
    reference_hirshfeld_weights_host(
      shifted_mol, shifted_meta, shifted_tasks.begin(), shifted_tasks.end()
    );

    for( size_t ipt = 0; ipt < points.size(); ++ipt ) {
      CHECK( shifted_tasks[0].weights[ipt] == Approx(tasks[0].weights[ipt]) );
    }
  }

  SECTION("atom order permutation preserves physical atom partitions") {
    Molecule mol;
    mol.emplace_back( AtomicNumber(6), -0.8,  0.0,  0.0 );
    mol.emplace_back( AtomicNumber(8),  1.3,  0.2, -0.1 );
    mol.emplace_back( AtomicNumber(1),  0.0, -1.1,  0.7 );

    Molecule permuted;
    permuted.push_back(mol[2]);
    permuted.push_back(mol[0]);
    permuted.push_back(mol[1]);
    const std::vector<size_t> original_to_permuted = {1, 2, 0};

    const std::vector<std::array<double,3>> points = {{
      {-0.2, -0.1,  0.1},
      { 1.0,  0.4, -0.3},
      { 0.2, -0.9,  0.6}
    }};

    MolMeta meta(mol);
    MolMeta permuted_meta(permuted);

    std::vector<XCTask> tasks(mol.natoms()), permuted_tasks(mol.natoms());
    for( size_t parent = 0; parent < mol.natoms(); ++parent ) {
      tasks[parent].iParent = static_cast<int32_t>(parent);
      tasks[parent].points = points;
      tasks[parent].weights = {1.0, 1.0, 1.0};
    }
    for( size_t parent = 0; parent < permuted.natoms(); ++parent ) {
      permuted_tasks[parent].iParent = static_cast<int32_t>(parent);
      permuted_tasks[parent].points = points;
      permuted_tasks[parent].weights = {1.0, 1.0, 1.0};
    }

    reference_hirshfeld_weights_host( mol, meta, tasks.begin(), tasks.end() );
    reference_hirshfeld_weights_host(
      permuted, permuted_meta, permuted_tasks.begin(), permuted_tasks.end()
    );

    for( size_t original_parent = 0; original_parent < mol.natoms(); ++original_parent ) {
      const auto permuted_parent = original_to_permuted[original_parent];
      for( size_t ipt = 0; ipt < points.size(); ++ipt ) {
        CHECK( permuted_tasks[permuted_parent].weights[ipt] ==
               Approx(tasks[original_parent].weights[ipt]) );
      }
    }
  }
}
#endif


