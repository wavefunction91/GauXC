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
#include "../src/xc_integrator/replicated/host/vv10_nlc.hpp"
#include "standards.hpp"

#include <gauxc/exceptions.hpp>
#include <gauxc/external/hdf5.hpp>
#include <gauxc/molgrid/defaults.hpp>
#include <gauxc/molecular_weights.hpp>
#include <gauxc/reduction_driver.hpp>
#include <gauxc/runtime_environment.hpp>
#include <gauxc/xc_integrator.hpp>
#include <gauxc/xc_integrator/impl.hpp>
#include <gauxc/xc_integrator/integrator_factory.hpp>

#include <Eigen/Core>
#include <highfive/H5File.hpp>

#include <array>
#include <cmath>
#include <numeric>

namespace {

std::vector<double> read_double_dataset(HighFive::File& file, const std::string& name) {
  auto dataset = file.getDataSet(name);
  auto dims = dataset.getDimensions();
  size_t size = 1;
  for( auto dim : dims ) size *= dim;
  std::vector<double> values(size);
  if( dims.size() == 1 ) {
    dataset.read(values);
  } else if( dims.size() == 2 ) {
    std::vector<std::vector<double>> matrix;
    dataset.read(matrix);
    size_t index = 0;
    for( const auto& row : matrix ) {
      for( const auto value : row ) values[index++] = value;
    }
  } else {
    FAIL("Unsupported VV10 reference dataset rank");
  }
  return values;
}

std::vector<unsigned long long> read_uint64_dataset(HighFive::File& file, const std::string& name) {
  auto dataset = file.getDataSet(name);
  auto dims = dataset.getDimensions();
  size_t size = 1;
  for( auto dim : dims ) size *= dim;
  std::vector<unsigned long long> values(size);
  dataset.read(values);
  return values;
}

}

TEST_CASE( "VV10 NLC helper matches PySCF reference", "[vv10][nlc]" ) {
  GauXC::IntegratorSettingsNLC settings;
  settings.vv10_b = 6.3;
  settings.vv10_c = 0.0093;

  constexpr std::size_t npts = 4;
  const std::array<double, 3*npts> coords = {
     0.00,  0.00,  0.00,
     0.37, -0.22,  0.51,
    -0.44,  0.18, -0.31,
     0.11,  0.42, -0.26
  };
  const std::array<double, npts> weights = { 0.31, 0.27, 0.23, 0.19 };
  const std::array<double, npts> rho = { 0.42, 0.37, 0.28, 0.51 };
  const std::array<double, npts> gamma = {
    0.030*0.030 + (-0.018)*(-0.018) + 0.022*0.022,
    (-0.020)*(-0.020) + 0.011*0.011 + 0.017*0.017,
    0.015*0.015 + 0.019*0.019 + (-0.012)*(-0.012),
    0.025*0.025 + (-0.014)*(-0.014) + 0.020*0.020
  };

  std::array<double, npts> eps = { 0.0, 0.0, 0.0, 0.0 };
  std::array<double, npts> vrho = { 0.0, 0.0, 0.0, 0.0 };
  std::array<double, npts> vgamma = { 0.0, 0.0, 0.0, 0.0 };

  GauXC::detail::vv10::eval_exc_vxc(
    settings,
    GauXC::detail::vv10::GridView{ npts, coords.data(), weights.data(), rho.data(), gamma.data() },
    GauXC::detail::vv10::CorrectionsView{ eps.data(), vrho.data(), vgamma.data() }
  );

  const std::array<double, npts> ref_eps = {
    0.0044608986634366,
    0.00446368137730377,
    0.00445955256959738,
    0.00446493156065464
  };
  const std::array<double, npts> ref_vrho = {
    0.00443982739412034,
    0.00444498672971984,
    0.00443833501920801,
    0.00444652164147545
  };
  const std::array<double, npts> ref_vgamma = {
    3.1648948742458134e-10,
    5.5585365145093665e-10,
    1.3162234665769539e-09,
    1.6830142122261283e-10
  };

  for( std::size_t i = 0; i < npts; ++i ) {
    CHECK( eps[i] == Approx(ref_eps[i]).margin(1e-14) );
    CHECK( vrho[i] == Approx(ref_vrho[i]).margin(1e-14) );
    CHECK( vgamma[i] == Approx(ref_vgamma[i]).margin(1e-18) );
  }
}

TEST_CASE( "VV10 NLC helper matches PySCF HDF5 reference fixture", "[vv10][nlc]" ) {
  HighFive::File file( GAUXC_REF_DATA_PATH "/vv10_pyscf_reference.hdf5", HighFive::File::ReadOnly );

  const auto params = read_double_dataset(file, "params");
  const auto coords = read_double_dataset(file, "coords");
  const auto weights = read_double_dataset(file, "weights");
  const auto rho = read_double_dataset(file, "rho");
  const auto grad = read_double_dataset(file, "grad");
  const auto gamma = read_double_dataset(file, "gamma");
  const auto rho_t = read_double_dataset(file, "rho_t");
  const auto grad_t = read_double_dataset(file, "grad_t");
  const auto parents = read_uint64_dataset(file, "parents");
  const auto ref_eps = read_double_dataset(file, "eps");
  const auto ref_vrho = read_double_dataset(file, "vrho");
  const auto ref_vgamma = read_double_dataset(file, "vgamma");
  const auto ref_response_A = read_double_dataset(file, "response_A");
  const auto ref_response_B = read_double_dataset(file, "response_B");
  const auto ref_grid_gradient = read_double_dataset(file, "grid_gradient");
  const auto ref_grid_gradient_excluding_same_parent = read_double_dataset(file, "grid_gradient_excluding_same_parent");

  const auto npts = rho.size();
  REQUIRE( params.size() == 2 );
  REQUIRE( coords.size() == 3 * npts );
  REQUIRE( weights.size() == npts );
  REQUIRE( grad.size() == 3 * npts );
  REQUIRE( grad_t.size() == 3 * npts );

  GauXC::IntegratorSettingsNLC settings;
  settings.vv10_b = params[0];
  settings.vv10_c = params[1];

  std::vector<double> eps(npts, 0.0), vrho(npts, 0.0), vgamma(npts, 0.0);
  GauXC::detail::vv10::eval_exc_vxc(
    settings,
    GauXC::detail::vv10::GridView{ npts, coords.data(), weights.data(), rho.data(), gamma.data() },
    GauXC::detail::vv10::CorrectionsView{ eps.data(), vrho.data(), vgamma.data() }
  );

  std::vector<double> grad_x(npts), grad_y(npts), grad_z(npts);
  std::vector<double> grad_t_x(npts), grad_t_y(npts), grad_t_z(npts);
  for( size_t i = 0; i < npts; ++i ) {
    grad_x[i] = grad[3*i];
    grad_y[i] = grad[3*i+1];
    grad_z[i] = grad[3*i+2];
    grad_t_x[i] = grad_t[3*i];
    grad_t_y[i] = grad_t[3*i+1];
    grad_t_z[i] = grad_t[3*i+2];
  }

  std::vector<double> response_A(npts, 0.0), response_B(3*npts, 0.0);
  GauXC::detail::vv10::eval_fxc_rks(
    settings,
    GauXC::detail::vv10::ResponseGridView{
      npts, coords.data(), weights.data(), rho.data(), gamma.data(), grad_x.data(), grad_y.data(), grad_z.data()
    },
    GauXC::detail::vv10::TrialView{ rho_t.data(), grad_t_x.data(), grad_t_y.data(), grad_t_z.data() },
    GauXC::detail::vv10::ResponseCorrectionsView{ response_A.data(), response_B.data() }
  );

  std::vector<double> grid_grad_x(npts, 0.0), grid_grad_y(npts, 0.0), grid_grad_z(npts, 0.0);
  GauXC::detail::vv10::eval_grid_gradient(
    settings,
    GauXC::detail::vv10::GridView{ npts, coords.data(), weights.data(), rho.data(), gamma.data() },
    GauXC::detail::vv10::GridGradientView{ grid_grad_x.data(), grid_grad_y.data(), grid_grad_z.data() }
  );

  std::vector<double> excluded_grad_x(npts, 0.0), excluded_grad_y(npts, 0.0), excluded_grad_z(npts, 0.0);
  GauXC::detail::vv10::eval_grid_gradient_excluding_same_parent(
    settings,
    GauXC::detail::vv10::GridView{ npts, coords.data(), weights.data(), rho.data(), gamma.data() },
    parents,
    GauXC::detail::vv10::GridGradientView{ excluded_grad_x.data(), excluded_grad_y.data(), excluded_grad_z.data() }
  );

  for( size_t i = 0; i < npts; ++i ) {
    CHECK( eps[i] == Approx(ref_eps[i]).margin(1e-14) );
    CHECK( vrho[i] == Approx(ref_vrho[i]).margin(1e-14) );
    CHECK( vgamma[i] == Approx(ref_vgamma[i]).margin(1e-18) );
    CHECK( response_A[i] == Approx(ref_response_A[i]).margin(1e-18) );
    CHECK( grid_grad_x[i] == Approx(ref_grid_gradient[3*i]).margin(1e-18) );
    CHECK( grid_grad_y[i] == Approx(ref_grid_gradient[3*i+1]).margin(1e-18) );
    CHECK( grid_grad_z[i] == Approx(ref_grid_gradient[3*i+2]).margin(1e-18) );
    CHECK( excluded_grad_x[i] == Approx(ref_grid_gradient_excluding_same_parent[3*i]).margin(1e-18) );
    CHECK( excluded_grad_y[i] == Approx(ref_grid_gradient_excluding_same_parent[3*i+1]).margin(1e-18) );
    CHECK( excluded_grad_z[i] == Approx(ref_grid_gradient_excluding_same_parent[3*i+2]).margin(1e-18) );
  }
  for( size_t i = 0; i < 3*npts; ++i ) {
    CHECK( response_B[i] == Approx(ref_response_B[i]).margin(1e-21) );
  }
}

TEST_CASE( "VV10 NLC UKS channel mapping", "[vv10][nlc]" ) {
  std::array<double, 2> eps = { 1.0, 2.0 };
  std::array<double, 4> vrho = { 0.0, 0.0, 10.0, 20.0 };
  std::array<double, 6> vgamma = { 0.0, 0.0, 0.0, 30.0, 40.0, 50.0 };

  GauXC::detail::vv10::add_uks_correction( 1, 0.25, 0.5, 0.75,
    eps.data(), vrho.data(), vgamma.data() );

  CHECK( eps[0] == Approx(1.0) );
  CHECK( eps[1] == Approx(2.25) );
  CHECK( vrho[0] == Approx(0.0) );
  CHECK( vrho[1] == Approx(0.0) );
  CHECK( vrho[2] == Approx(10.5) );
  CHECK( vrho[3] == Approx(20.5) );
  CHECK( vgamma[0] == Approx(0.0) );
  CHECK( vgamma[1] == Approx(0.0) );
  CHECK( vgamma[2] == Approx(0.0) );
  CHECK( vgamma[3] == Approx(30.75) );
  CHECK( vgamma[4] == Approx(41.5) );
  CHECK( vgamma[5] == Approx(50.75) );
}

TEST_CASE( "VV10 NLC RKS response matches PySCF libdft reference", "[vv10][nlc]" ) {
  GauXC::IntegratorSettingsNLC settings;
  settings.vv10_b = 6.3;
  settings.vv10_c = 0.0093;

  constexpr std::size_t npts = 4;
  const std::array<double, 3*npts> coords = {
     0.00,  0.00,  0.00,
     0.37, -0.22,  0.51,
    -0.44,  0.18, -0.31,
     0.11,  0.42, -0.26
  };
  const std::array<double, npts> weights = { 0.31, 0.27, 0.23, 0.19 };
  const std::array<double, npts> rho = { 0.42, 0.37, 0.28, 0.51 };
  const std::array<double, npts> grad_x = { 0.030, -0.020, 0.015, 0.025 };
  const std::array<double, npts> grad_y = { -0.018, 0.011, 0.019, -0.014 };
  const std::array<double, npts> grad_z = { 0.022, 0.017, -0.012, 0.020 };
  const std::array<double, npts> gamma = {
    0.030*0.030 + (-0.018)*(-0.018) + 0.022*0.022,
    (-0.020)*(-0.020) + 0.011*0.011 + 0.017*0.017,
    0.015*0.015 + 0.019*0.019 + (-0.012)*(-0.012),
    0.025*0.025 + (-0.014)*(-0.014) + 0.020*0.020
  };
  const std::array<double, npts> rho_t = { 0.013, -0.021, 0.017, 0.009 };
  const std::array<double, npts> grad_t_x = { 0.004, -0.003, 0.002, 0.006 };
  const std::array<double, npts> grad_t_y = { -0.002, 0.005, -0.004, 0.001 };
  const std::array<double, npts> grad_t_z = { 0.003, 0.002, 0.005, -0.003 };

  std::array<double, npts> A = { 0.0, 0.0, 0.0, 0.0 };
  std::array<double, 3*npts> B = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  GauXC::detail::vv10::eval_fxc_rks(
    settings,
    GauXC::detail::vv10::ResponseGridView{
      npts, coords.data(), weights.data(), rho.data(), gamma.data(), grad_x.data(), grad_y.data(), grad_z.data()
    },
    GauXC::detail::vv10::TrialView{ rho_t.data(), grad_t_x.data(), grad_t_y.data(), grad_t_z.data() },
    GauXC::detail::vv10::ResponseCorrectionsView{ A.data(), B.data() }
  );

  const std::array<double, npts> ref_A = {
    -3.5703422236383133e-08,
    -1.1793774790700641e-06,
     2.9325704002469295e-07,
    -3.5681576663273375e-07
  };
  const std::array<double, 3*npts> ref_B = {
     5.0726028793207675e-12,
    -2.7903701376527950e-12,
     3.7621073764918405e-12,
    -1.7056375465935961e-11,
     1.3105225970986051e-11,
     1.3886480129449537e-11,
    -1.6237863203856369e-11,
    -3.7766613354823595e-11,
     3.0364440321900880e-11,
     2.3991229362623262e-12,
     1.2407954875428006e-13,
    -7.0620382206289707e-13
  };

  for( std::size_t i = 0; i < npts; ++i ) {
    CHECK( A[i] == Approx(ref_A[i]).margin(1e-18) );
  }
  for( std::size_t i = 0; i < 3*npts; ++i ) {
    CHECK( B[i] == Approx(ref_B[i]).margin(1e-21) );
  }
}

TEST_CASE( "VV10 NLC grid gradient matches PySCF reference", "[vv10][nlc]" ) {
  GauXC::IntegratorSettingsNLC settings;
  settings.vv10_b = 6.3;
  settings.vv10_c = 0.0093;

  constexpr std::size_t npts = 4;
  const std::array<double, 3*npts> coords = {
     0.00,  0.00,  0.00,
     0.37, -0.22,  0.51,
    -0.44,  0.18, -0.31,
     0.11,  0.42, -0.26
  };
  const std::array<double, npts> weights = { 0.31, 0.27, 0.23, 0.19 };
  const std::array<double, npts> rho = { 0.42, 0.37, 0.28, 0.51 };
  const std::array<double, npts> gamma = {
    0.030*0.030 + (-0.018)*(-0.018) + 0.022*0.022,
    (-0.020)*(-0.020) + 0.011*0.011 + 0.017*0.017,
    0.015*0.015 + 0.019*0.019 + (-0.012)*(-0.012),
    0.025*0.025 + (-0.014)*(-0.014) + 0.020*0.020
  };
  std::array<double, npts> grad_x = { 0.0, 0.0, 0.0, 0.0 };
  std::array<double, npts> grad_y = { 0.0, 0.0, 0.0, 0.0 };
  std::array<double, npts> grad_z = { 0.0, 0.0, 0.0, 0.0 };

  GauXC::detail::vv10::eval_grid_gradient(
    settings,
    GauXC::detail::vv10::GridView{ npts, coords.data(), weights.data(), rho.data(), gamma.data() },
    GauXC::detail::vv10::GridGradientView{ grad_x.data(), grad_y.data(), grad_z.data() }
  );

  const std::array<double, 3*npts> ref_grad = {
    -1.9307886424617078e-06,
    -3.5876026873695561e-06,
    -3.7865972833005390e-07,
     1.1869243697781221e-05,
    -1.0719376926550508e-05,
     1.8240846117598208e-05,
    -1.9673087097543844e-05,
     3.4494080789274471e-06,
    -1.2190661056077449e-05,
     3.4323843645200651e-06,
     1.3579254329979196e-05,
    -1.0194844773045410e-05
  };

  for( std::size_t i = 0; i < npts; ++i ) {
    CHECK( grad_x[i] == Approx(ref_grad[3*i]).margin(1e-18) );
    CHECK( grad_y[i] == Approx(ref_grad[3*i+1]).margin(1e-18) );
    CHECK( grad_z[i] == Approx(ref_grad[3*i+2]).margin(1e-18) );
  }
}

TEST_CASE( "VV10 NLC grid gradient can exclude same-parent grid points", "[vv10][nlc]" ) {
  GauXC::IntegratorSettingsNLC settings;
  settings.vv10_b = 6.3;
  settings.vv10_c = 0.0093;

  constexpr std::size_t npts = 4;
  const std::array<double, 3*npts> coords = {
     0.00,  0.00,  0.00,
     0.37, -0.22,  0.51,
    -0.44,  0.18, -0.31,
     0.11,  0.42, -0.26
  };
  const std::array<double, npts> weights = { 0.31, 0.27, 0.23, 0.19 };
  const std::array<double, npts> rho = { 0.42, 0.37, 0.28, 0.51 };
  const std::array<double, npts> gamma = {
    0.030*0.030 + (-0.018)*(-0.018) + 0.022*0.022,
    (-0.020)*(-0.020) + 0.011*0.011 + 0.017*0.017,
    0.015*0.015 + 0.019*0.019 + (-0.012)*(-0.012),
    0.025*0.025 + (-0.014)*(-0.014) + 0.020*0.020
  };
  std::array<double, npts> grad_x = { 0.0, 0.0, 0.0, 0.0 };
  std::array<double, npts> grad_y = { 0.0, 0.0, 0.0, 0.0 };
  std::array<double, npts> grad_z = { 0.0, 0.0, 0.0, 0.0 };
  const std::vector<unsigned long long> parent = { 0, 1, 1, 0 };

  GauXC::detail::vv10::eval_grid_gradient_excluding_same_parent(
    settings,
    GauXC::detail::vv10::GridView{ npts, coords.data(), weights.data(), rho.data(), gamma.data() },
    parent,
    GauXC::detail::vv10::GridGradientView{ grad_x.data(), grad_y.data(), grad_z.data() }
  );

  const std::array<double, 3*npts> ref_grad = {
    -7.1657342952056772e-07,
     1.0484917620420706e-06,
    -3.2486229589182034e-06,
     7.5039909408315638e-06,
    -8.5636965527482088e-06,
     1.3821701351303496e-05,
    -1.2901522650039642e-05,
     1.0542563571549509e-07,
    -5.3354970474929469e-06,
     1.8009001465124650e-06,
     7.3499509521319973e-06,
    -6.3386093486638098e-06
  };

  for( std::size_t i = 0; i < npts; ++i ) {
    CHECK( grad_x[i] == Approx(ref_grad[3*i]).margin(1e-18) );
    CHECK( grad_y[i] == Approx(ref_grad[3*i+1]).margin(1e-18) );
    CHECK( grad_z[i] == Approx(ref_grad[3*i+2]).margin(1e-18) );
  }
}

TEST_CASE( "VV10 packed-grid allgather handles variable local sizes", "[vv10][nlc]" ) {
  auto rt = GauXC::RuntimeEnvironment(GAUXC_MPI_CODE(MPI_COMM_WORLD));
  auto reduction_driver = GauXC::ReductionDriverFactory::get_shared_instance(rt, "Default");

  const auto rank = reduction_driver->comm_rank();
  const auto size = reduction_driver->comm_size();
  std::vector<double> local( 2 * static_cast<std::size_t>(rank + 1) );
  for( std::size_t i = 0; i < local.size(); ++i ) {
    local[i] = 100.0 * rank + static_cast<double>(i);
  }

  std::size_t local_point_offset = 0;
  const auto gathered = GauXC::detail::vv10::allgather_packed_grid(
    *reduction_driver, local, 2, local_point_offset );

  std::size_t expected_points = 0;
  for( int i = 0; i < size; ++i ) expected_points += static_cast<std::size_t>(i + 1);
  CHECK( gathered.size() == 2 * expected_points );

  std::size_t expected_offset = 0;
  for( int i = 0; i < rank; ++i ) expected_offset += static_cast<std::size_t>(i + 1);
  CHECK( local_point_offset == expected_offset );

  std::size_t point_offset = 0;
  for( int irank = 0; irank < size; ++irank ) {
    for( int ipt = 0; ipt < irank + 1; ++ipt ) {
      CHECK( gathered[2*(point_offset + ipt)] == Approx(100.0 * irank + 2.0 * ipt) );
      CHECK( gathered[2*(point_offset + ipt)+1] == Approx(100.0 * irank + 2.0 * ipt + 1.0) );
    }
    point_offset += static_cast<std::size_t>(irank + 1);
  }
}

TEST_CASE( "VV10 integrated host paths are finite and repeatable", "[vv10][integrated]" ) {
  using matrix_type = Eigen::MatrixXd;

  auto rt = GauXC::RuntimeEnvironment(GAUXC_MPI_CODE(MPI_COMM_WORLD));
  auto mol = GauXC::make_water();
  auto basis = GauXC::make_631Gd( mol, GauXC::SphericalType(true) );
  for( auto& shell : basis ) {
    shell.set_shell_tolerance( std::numeric_limits<double>::epsilon() );
  }

  auto molgrid = GauXC::MolGridFactory::create_default_molgrid(
    mol, GauXC::PruningScheme::Unpruned, GauXC::BatchSize(128),
    GauXC::RadialQuad::MuraKnowles, GauXC::AtomicGridSizeDefault::GM3 );

  GauXC::LoadBalancerFactory lb_factory( GauXC::ExecutionSpace::Host, "Default" );
  auto load_balancer = lb_factory.get_instance( rt, mol, molgrid, basis );

  GauXC::MolecularWeightsFactory mw_factory(
    GauXC::ExecutionSpace::Host, "Default", GauXC::MolecularWeightsSettings{} );
  auto molecular_weights = mw_factory.get_instance();
  molecular_weights.modify_weights( load_balancer );

  auto functional = GauXC::functional_type(
    ExchCXX::Backend::builtin, ExchCXX::Functional::BLYP, ExchCXX::Spin::Unpolarized );
  GauXC::XCIntegratorFactory<matrix_type> integrator_factory(
    GauXC::ExecutionSpace::Host, "Replicated", "Reference", "Default", "Default" );
  auto integrator = integrator_factory.get_instance( functional, load_balancer );

  const auto nbf = basis.nbf();
  matrix_type P = matrix_type::Zero( nbf, nbf );
  matrix_type tP = matrix_type::Zero( nbf, nbf );
  for( int i = 0; i < nbf; ++i ) {
    P(i,i) = 0.05 + 0.002 * i;
    tP(i,i) = 0.01 + 0.001 * i;
  }

  GauXC::IntegratorSettingsNLC settings;
  settings.include_weight_derivatives = true;

  auto [base_exc, base_vxc] = integrator.eval_exc_vxc( P );
  auto base_fxc = integrator.eval_fxc_contraction( P, tP );

  auto [exc, vxc] = integrator.eval_exc_vxc( P, settings );
  auto [exc_repeat, vxc_repeat] = integrator.eval_exc_vxc( P, settings );
  CHECK( std::isfinite(exc) );
  CHECK( exc == Approx(base_exc).margin(1e-12) );
  CHECK( (vxc - base_vxc).norm() < 1e-10 );
  CHECK( exc_repeat == Approx(base_exc).margin(1e-12) );
  CHECK( (vxc - vxc_repeat).norm() < 1e-10 );
  CHECK( integrator.eval_exc( P, settings ) == Approx(base_exc).margin(1e-10) );

  auto [enlc, vnlc] = integrator.eval_nlc_vnlc( P, settings );
  CHECK( std::isfinite(enlc) );
  CHECK( std::abs(enlc) > 1e-10 );
  CHECK( vnlc.allFinite() );
  CHECK( vnlc.norm() > 1e-10 );
  CHECK( integrator.eval_nlc( P, settings ) == Approx(enlc).margin(1e-10) );

  auto exc_grad_full = integrator.eval_exc_grad( P, settings );
  auto nlc_exc_grad_full = integrator.eval_nlc_grad( P, settings );
  settings.include_weight_derivatives = false;
  auto exc_grad_hf = integrator.eval_exc_grad( P, settings );
  auto nlc_exc_grad_hf = integrator.eval_nlc_grad( P, settings );
  REQUIRE( exc_grad_full.size() == 3 * mol.natoms() );
  REQUIRE( exc_grad_hf.size() == 3 * mol.natoms() );
  REQUIRE( nlc_exc_grad_full.size() == 3 * mol.natoms() );
  REQUIRE( nlc_exc_grad_hf.size() == 3 * mol.natoms() );
  for( auto value : exc_grad_full ) CHECK( std::isfinite(value) );
  for( auto value : exc_grad_hf ) CHECK( std::isfinite(value) );
  for( auto value : nlc_exc_grad_full ) CHECK( std::isfinite(value) );
  for( auto value : nlc_exc_grad_hf ) CHECK( std::isfinite(value) );
  CHECK( std::sqrt(std::inner_product(nlc_exc_grad_full.begin(), nlc_exc_grad_full.end(), nlc_exc_grad_full.begin(), 0.0)) > 1e-12 );
  CHECK( std::sqrt(std::inner_product(nlc_exc_grad_hf.begin(), nlc_exc_grad_hf.end(), nlc_exc_grad_hf.begin(), 0.0)) > 1e-12 );

  double full_sum[3] = { 0.0, 0.0, 0.0 };
  for( size_t iatom = 0; iatom < mol.natoms(); ++iatom ) {
    full_sum[0] += exc_grad_full[3*iatom+0];
    full_sum[1] += exc_grad_full[3*iatom+1];
    full_sum[2] += exc_grad_full[3*iatom+2];
  }
  CHECK( std::abs(full_sum[0]) < 1e-8 );
  CHECK( std::abs(full_sum[1]) < 1e-8 );
  CHECK( std::abs(full_sum[2]) < 1e-8 );

  double nlc_full_sum[3] = { 0.0, 0.0, 0.0 };
  for( size_t iatom = 0; iatom < mol.natoms(); ++iatom ) {
    nlc_full_sum[0] += nlc_exc_grad_full[3*iatom+0];
    nlc_full_sum[1] += nlc_exc_grad_full[3*iatom+1];
    nlc_full_sum[2] += nlc_exc_grad_full[3*iatom+2];
  }
  CHECK( std::abs(nlc_full_sum[0]) < 1e-8 );
  CHECK( std::abs(nlc_full_sum[1]) < 1e-8 );
  CHECK( std::abs(nlc_full_sum[2]) < 1e-8 );

  auto fxc = integrator.eval_fxc_contraction( P, tP, settings );
  CHECK( fxc.allFinite() );
  CHECK( (fxc - fxc.transpose()).norm() / nbf < 1e-10 );
  CHECK( (fxc - base_fxc).norm() < 1e-10 );
  auto nlc_fxc = integrator.eval_nlc_fnlc_contraction( P, tP, settings );
  CHECK( nlc_fxc.allFinite() );
  CHECK( (nlc_fxc - nlc_fxc.transpose()).norm() / nbf < 1e-10 );
  CHECK( nlc_fxc.norm() > 1e-12 );

  matrix_type Pz = matrix_type::Zero( nbf, nbf );
  matrix_type tPz = matrix_type::Zero( nbf, nbf );
  auto polarized_functional = GauXC::functional_type(
    ExchCXX::Backend::builtin, ExchCXX::Functional::BLYP, ExchCXX::Spin::Polarized );
  auto uks_integrator = integrator_factory.get_instance( polarized_functional, load_balancer );
  auto [uks_base_exc, uks_base_vxc, uks_base_vxcz] = uks_integrator.eval_exc_vxc( P, Pz );
  auto [uks_exc, uks_vxc, uks_vxcz] = uks_integrator.eval_exc_vxc( P, Pz, settings );
  auto [uks_enlc, uks_vnlc, uks_vnlcz] = uks_integrator.eval_nlc_vnlc( P, Pz, settings );
  CHECK( std::isfinite(uks_exc) );
  CHECK( uks_exc == Approx(uks_base_exc).margin(1e-12) );
  CHECK( (uks_vxc - uks_base_vxc).norm() < 1e-10 );
  CHECK( (uks_vxcz - uks_base_vxcz).norm() < 1e-10 );
  CHECK( std::isfinite(uks_enlc) );
  CHECK( std::abs(uks_enlc) > 1e-10 );
  // A spin-unpolarized UKS density must reproduce the RKS NLC energy exactly,
  // since VV10 depends only on the total density and its gradient. RKS scales
  // the density matrix by 2 (closed-shell double occupancy) while UKS treats
  // the scalar matrix as the full total density, so the equivalent UKS input
  // is 2*P with zero spin density.
  matrix_type P2 = 2.0 * P;
  CHECK( uks_integrator.eval_nlc( P2, Pz, settings ) == Approx(enlc).margin(1e-10) );
  CHECK( uks_vnlc.allFinite() );
  CHECK( uks_vnlc.norm() > 1e-10 );
  CHECK( uks_vnlcz.norm() / nbf < 1e-8 );
  CHECK( uks_integrator.eval_nlc( P, Pz, settings ) == Approx(uks_enlc).margin(1e-10) );
  CHECK( uks_vxc.allFinite() );
  CHECK( uks_vxcz.allFinite() );
  CHECK( uks_vxcz.norm() / nbf < 1e-8 );
  auto uks_nlc_grad = uks_integrator.eval_nlc_grad( P, Pz, settings );
  REQUIRE( uks_nlc_grad.size() == 3 * mol.natoms() );
  for( auto value : uks_nlc_grad ) CHECK( std::isfinite(value) );
  CHECK( std::sqrt(std::inner_product(uks_nlc_grad.begin(), uks_nlc_grad.end(), uks_nlc_grad.begin(), 0.0)) > 1e-12 );
  auto [uks_base_fxc, uks_base_fxcz] = uks_integrator.eval_fxc_contraction( P, Pz, tP, tPz );
  auto [uks_fxc, uks_fxcz] = uks_integrator.eval_fxc_contraction( P, Pz, tP, tPz, settings );
  auto [uks_nlc_fxc, uks_nlc_fxcz] = uks_integrator.eval_nlc_fnlc_contraction( P, Pz, tP, tPz, settings );
  CHECK( uks_fxc.allFinite() );
  CHECK( uks_fxcz.allFinite() );
  CHECK( (uks_fxc - uks_base_fxc).norm() < 1e-10 );
  CHECK( (uks_fxcz - uks_base_fxcz).norm() < 1e-10 );
  CHECK( uks_fxcz.norm() / nbf < 1e-8 );
  CHECK( uks_nlc_fxc.allFinite() );
  CHECK( uks_nlc_fxc.norm() > 1e-12 );
  CHECK( uks_nlc_fxcz.norm() / nbf < 1e-8 );
}

#ifdef GAUXC_HAS_DEVICE
TEST_CASE( "VV10 integrated device RKS EXC/VXC matches host", "[vv10][integrated][cuda]" ) {
  using matrix_type = Eigen::MatrixXd;

  auto host_rt = GauXC::RuntimeEnvironment(GAUXC_MPI_CODE(MPI_COMM_WORLD));
  auto device_rt = GauXC::DeviceRuntimeEnvironment(GAUXC_MPI_CODE(MPI_COMM_WORLD,) 0.5);
  auto mol = GauXC::make_water();
  auto basis = GauXC::make_631Gd( mol, GauXC::SphericalType(true) );
  for( auto& shell : basis ) {
    shell.set_shell_tolerance( std::numeric_limits<double>::epsilon() );
  }

  auto molgrid = GauXC::MolGridFactory::create_default_molgrid(
    mol, GauXC::PruningScheme::Unpruned, GauXC::BatchSize(128),
    GauXC::RadialQuad::MuraKnowles, GauXC::AtomicGridSizeDefault::GM3 );

  GauXC::LoadBalancerFactory lb_factory( GauXC::ExecutionSpace::Host, "Default" );
  auto host_lb = lb_factory.get_instance( host_rt, mol, molgrid, basis );
  auto device_lb = lb_factory.get_instance( device_rt, mol, molgrid, basis );

  GauXC::MolecularWeightsFactory host_mw_factory(
    GauXC::ExecutionSpace::Host, "Default", GauXC::MolecularWeightsSettings{} );
  auto host_mw = host_mw_factory.get_instance();
  host_mw.modify_weights( host_lb );

  GauXC::MolecularWeightsFactory device_mw_factory(
    GauXC::ExecutionSpace::Device, "Default", GauXC::MolecularWeightsSettings{} );
  auto device_mw = device_mw_factory.get_instance();
  device_mw.modify_weights( device_lb );

  auto functional = GauXC::functional_type(
    ExchCXX::Backend::builtin, ExchCXX::Functional::BLYP, ExchCXX::Spin::Unpolarized );
  GauXC::XCIntegratorFactory<matrix_type> host_integrator_factory(
    GauXC::ExecutionSpace::Host, "Replicated", "Reference", "Default", "Default" );
  auto host_integrator = host_integrator_factory.get_instance( functional, host_lb );
  GauXC::XCIntegratorFactory<matrix_type> device_integrator_factory(
    GauXC::ExecutionSpace::Device, "Replicated", "Default", "Default", "Default" );
  auto device_integrator = device_integrator_factory.get_instance( functional, device_lb );

  const auto nbf = basis.nbf();
  matrix_type P = matrix_type::Zero( nbf, nbf );
  matrix_type tP = matrix_type::Zero( nbf, nbf );
  for( int i = 0; i < nbf; ++i ) {
    P(i,i) = 0.05 + 0.002 * i;
    tP(i,i) = 0.01 + 0.001 * i;
  }

  GauXC::IntegratorSettingsNLC settings;
  auto [host_enlc, host_vnlc] = host_integrator.eval_nlc_vnlc( P, settings );
  auto [device_enlc, device_vnlc] = device_integrator.eval_nlc_vnlc( P, settings );

  CHECK( device_enlc == Approx(host_enlc).margin(1e-8) );
  CHECK( (device_vnlc - host_vnlc).norm() / nbf < 1e-7 );
  CHECK( (device_vnlc - device_vnlc.transpose()).norm() / nbf < 1e-10 );
  CHECK( device_integrator.eval_nlc( P, settings ) == Approx(host_enlc).margin(1e-8) );
  matrix_type Pz = matrix_type::Zero( nbf, nbf );
  matrix_type tPz = matrix_type::Zero( nbf, nbf );
  auto polarized_functional = GauXC::functional_type(
    ExchCXX::Backend::builtin, ExchCXX::Functional::BLYP, ExchCXX::Spin::Polarized );
  auto polarized_host_integrator = host_integrator_factory.get_instance( polarized_functional, host_lb );
  auto polarized_device_integrator = device_integrator_factory.get_instance( polarized_functional, device_lb );
  auto [host_uks_enlc, host_uks_vnlc, host_uks_vnlcz] = polarized_host_integrator.eval_nlc_vnlc( P, Pz, settings );
  auto [device_uks_enlc, device_uks_vnlc, device_uks_vnlcz] = polarized_device_integrator.eval_nlc_vnlc( P, Pz, settings );
  CHECK( device_uks_enlc == Approx(host_uks_enlc).margin(1e-8) );
  CHECK( (device_uks_vnlc - host_uks_vnlc).norm() / nbf < 1e-7 );
  CHECK( device_uks_vnlcz.norm() / nbf < 1e-10 );
  CHECK( polarized_device_integrator.eval_nlc( P, Pz, settings ) == Approx(host_uks_enlc).margin(1e-8) );
  auto host_nlc_fxc = host_integrator.eval_nlc_fnlc_contraction( P, tP, settings );
  auto device_nlc_fxc = device_integrator.eval_nlc_fnlc_contraction( P, tP, settings );
  CHECK( device_nlc_fxc.allFinite() );
  CHECK( (device_nlc_fxc - host_nlc_fxc).norm() / nbf < 1e-7 );
  CHECK( (device_nlc_fxc - device_nlc_fxc.transpose()).norm() / nbf < 1e-10 );
  auto hf_settings = settings;
  hf_settings.include_weight_derivatives = false;
  auto host_nlc_grad = host_integrator.eval_nlc_grad( P, hf_settings );
  auto device_nlc_grad = device_integrator.eval_nlc_grad( P, hf_settings );
  REQUIRE( device_nlc_grad.size() == host_nlc_grad.size() );
  double grad_diff = 0.0;
  for( size_t i = 0; i < host_nlc_grad.size(); ++i ) {
    grad_diff += (device_nlc_grad[i] - host_nlc_grad[i]) * (device_nlc_grad[i] - host_nlc_grad[i]);
  }
  CHECK( std::sqrt(grad_diff) / mol.natoms() < 1e-7 );
  host_nlc_grad = host_integrator.eval_nlc_grad( P, settings );
  device_nlc_grad = device_integrator.eval_nlc_grad( P, settings );
  grad_diff = 0.0;
  for( size_t i = 0; i < host_nlc_grad.size(); ++i ) {
    grad_diff += (device_nlc_grad[i] - host_nlc_grad[i]) * (device_nlc_grad[i] - host_nlc_grad[i]);
  }
  CHECK( std::sqrt(grad_diff) / mol.natoms() < 1e-7 );
  CHECK_THROWS_AS( polarized_device_integrator.eval_nlc_fnlc_contraction( P, Pz, tP, tPz, settings ), GauXC::generic_gauxc_exception );
}
#endif