#include <catch2/catch.hpp>

#include <gauxc/gauxc_config.hpp>
#include <gauxc/reduction_driver.hpp>
#include <gauxc/runtime_environment.hpp>
#include "integrator_util/vv10.hpp"

#include <array>
#include <vector>

#ifdef GAUXC_HAS_CUDA
#include "device/cuda/kernels/vv10.hpp"
#include "device_specific/cuda_util.hpp"
#include "../src/xc_integrator/replicated/device/vv10_nlc_device.hpp"
#endif

namespace {

constexpr std::int64_t npts = 5;
constexpr std::int64_t ntrial = 2;
constexpr std::int64_t natoms = 2;

const std::vector<double> coords = {
   0.0,  0.1, -0.2,
   0.4, -0.3,  0.2,
  -0.5,  0.2,  0.3,
   0.7,  0.6, -0.4,
  -0.2, -0.7,  0.5
};

const std::vector<double> weights = {0.22, 0.17, 0.31, 0.11, 0.19};
const std::vector<double> rho = {0.72, 0.44, 0.91, 0.36, 0.58};
const std::vector<double> density_gradient_norm = {0.0089, 0.0078, 0.0158, 0.0101, 0.0125};
const std::vector<std::int32_t> grid_atoms = {0, 0, 1, 1, 0};
const std::vector<double> rho_t = {0.03, -0.02, 0.04, 0.01, -0.05, -0.01, 0.06, -0.03, 0.02, 0.04};
const std::vector<double> gamma_t = {0.002, -0.001, 0.003, -0.002, 0.001, -0.004, 0.002, 0.001, -0.003, 0.005};

template <typename T>
void check_vector(const std::vector<T>& actual, const std::vector<T>& ref, double margin = 1e-13) {
  REQUIRE(actual.size() == ref.size());
  for(std::size_t i = 0; i < actual.size(); ++i) {
    CHECK(actual[i] == Approx(ref[i]).margin(margin));
  }
}

std::vector<double> rho_weight() {
  std::vector<double> out(npts);
  for(std::int64_t i = 0; i < npts; ++i) out[i] = rho[i] * weights[i];
  return out;
}

}

TEST_CASE("VV10 host EXC/VXC formula matches PySCF reference", "[vv10][nlc]") {
  const GauXC::detail::vv10::Parameters params{6.3, 0.0093, 1e-8};
  std::vector<double> eps(npts), vrho(npts), vgamma(npts);

  GauXC::detail::vv10::eval_exc_vxc(npts, coords.data(), weights.data(),
    rho.data(), density_gradient_norm.data(), params, eps.data(), vrho.data(), vgamma.data());

  check_vector(eps, {
    0.00445301813026348, 0.00444920198423763, 0.00445616308451667,
    0.00445820144472496, 0.00445560260974276
  });
  check_vector(vrho, {
    0.00442962367367506, 0.00442417961353026, 0.00443444131349190,
    0.00443937740576438, 0.00443445791849522
  });
  check_vector(vgamma, {
    5.4703231981970949e-10, 4.3982859960426448e-09, 3.7195930885603543e-10,
    2.0026596148858508e-08, 2.6546741400669848e-09
  }, 1e-17);
}

TEST_CASE("VV10 variable-size host allgather support", "[vv10][nlc][mpi]") {
  const auto rt = GauXC::RuntimeEnvironment(GAUXC_MPI_CODE(MPI_COMM_WORLD));
  auto reduction = GauXC::ReductionDriverFactory::get_instance(rt, "BASICMPI");

  const auto rank = rt.comm_rank();
  const auto world_size = rt.comm_size();
  std::vector<double> local(rank + 1);
  for(int i = 0; i <= rank; ++i) local[i] = 100.0 * rank + i;

  auto gathered = reduction.allgather_v(local.data(), local.size());
  std::vector<double> reference;
  for(int r = 0; r < world_size; ++r) {
    for(int i = 0; i <= r; ++i) reference.push_back(100.0 * r + i);
  }

  check_vector(gathered, reference);
}

TEST_CASE("VV10 host gradient and response intermediates match references", "[vv10][nlc]") {
  const GauXC::detail::vv10::Parameters params{6.3, 0.0093, 1e-8};
  std::vector<double> omega(npts), kappa(npts), domega_drho(npts), domega_dgamma(npts);
  std::vector<double> d2omega_drho2(npts), d2omega_dgamma2(npts), d2omega_drho_dgamma(npts);
  std::vector<double> dkappa_drho(npts), d2kappa_drho2(npts);
  std::vector<double> U(npts), W(npts), A(npts), B(npts), C(npts), E(npts);

  GauXC::detail::vv10::eval_omega_kappa(npts, rho.data(), density_gradient_norm.data(), params,
    omega.data(), kappa.data(), domega_drho.data(), domega_dgamma.data(),
    d2omega_drho2.data(), d2omega_dgamma2.data(), d2omega_drho_dgamma.data(),
    dkappa_drho.data(), d2kappa_drho2.data());

  check_vector(omega, {1.7366438001493870, 1.3576018510928551, 1.9523837921688731, 1.2280150473242860, 1.5586889232762682});
  check_vector(kappa, {16.102934766228300, 14.834004175743555, 16.743895938069180, 14.346083941530297, 15.532962118860626});

  GauXC::detail::vv10::eval_hessian_intermediates(npts, coords.data(), weights.data(),
    rho.data(), omega.data(), kappa.data(), U.data(), W.data(), A.data(), B.data(), C.data(), E.data());

  check_vector(U, {9.209072244410352e-06, 1.046141581842596e-05, 8.359163476842639e-06, 8.499921019698521e-06, 8.763966041169824e-06});
  check_vector(W, {4.2839882446393177e-06, 7.0119105202783746e-06, 3.7243009447981262e-06, 1.2215577641798923e-05, 6.9448310913784774e-06});
  check_vector(E, {-1.0365967086190920e-04, -1.1129196291360637e-04, -9.7369762355538814e-05, -9.3293041938961706e-05, -9.8490711903360671e-05});

  std::vector<double> f_rho_t(ntrial*npts), f_gamma_t(ntrial*npts);
  GauXC::detail::vv10::eval_hessian_contraction(npts, ntrial, coords.data(), weights.data(),
    rho.data(), omega.data(), kappa.data(), U.data(), W.data(), A.data(), B.data(), C.data(),
    domega_drho.data(), domega_dgamma.data(), dkappa_drho.data(), d2omega_drho2.data(),
    d2omega_dgamma2.data(), d2omega_drho_dgamma.data(), d2kappa_drho2.data(), rho_t.data(),
    gamma_t.data(), f_rho_t.data(), f_gamma_t.data());

  check_vector(f_rho_t, {
     1.6178037606933446e-07, -1.3249866968720368e-06,  1.1116302821324279e-07,
    -7.5930158108261938e-08, -1.8068421774630219e-06, -1.0423414119308089e-06,
     1.8870904199228002e-06, -1.0109141134429277e-06,  3.5175453056695153e-07,
     3.5357534022654595e-07
  });
  check_vector(f_gamma_t, {
     2.8164183798520921e-11,  2.8193667944806545e-10,  2.2652562869246210e-13,
    -6.1629930379245268e-09,  1.1756835860583153e-09, -2.0337406941186164e-10,
    -1.2427741534208510e-09,  8.7638580360613370e-11, -1.0186539929376707e-08,
     3.4109976756300954e-10
  }, 1e-17);

  std::vector<double> Egr(natoms*3*npts), Ugr(natoms*3*npts), Wgr(natoms*3*npts);
  GauXC::detail::vv10::eval_grid_response(npts, natoms, coords.data(), weights.data(),
    rho.data(), omega.data(), kappa.data(), grid_atoms.data(), Egr.data(), Ugr.data(), Wgr.data());

  check_vector(Egr, {
     1.0310153822420043e-05,  1.8889529270293274e-05,  1.5112372095644556e-05,
    -1.6415873674929635e-05,  5.1457932465049436e-06, -4.4092456971783141e-06,
    -1.4045394092983390e-05, -1.2200427500773054e-05, -1.9148451277680992e-05,
    -2.2241372703956841e-05, -1.2143915488587469e-05, -2.1225565415921502e-07,
    -6.0778553859699155e-06,  1.0967407941688859e-05,  5.9820261926690516e-06,
    -1.0310153822420043e-05, -1.8889529270293274e-05, -1.5112372095644556e-05,
     1.6415873674929635e-05, -5.1457932465049436e-06,  4.4092456971783141e-06,
     1.4045394092983390e-05,  1.2200427500773054e-05,  1.9148451277680992e-05,
     2.2241372703956841e-05,  1.2143915488587469e-05,  2.1225565415921502e-07,
     6.0778553859699155e-06, -1.0967407941688859e-05, -5.9820261926690516e-06
  });
}

#ifdef GAUXC_HAS_CUDA
TEST_CASE("VV10 CUDA kernels match host formulas", "[vv10][nlc][cuda]") {
  const GauXC::detail::vv10::Parameters params{6.3, 0.0093, 1e-8};
  auto rw = rho_weight();

  auto* coords_d = GauXC::util::cuda_malloc<double>(coords.size());
  auto* weights_d = GauXC::util::cuda_malloc<double>(weights.size());
  auto* rho_d = GauXC::util::cuda_malloc<double>(rho.size());
  auto* gamma_d = GauXC::util::cuda_malloc<double>(density_gradient_norm.size());
  auto* rw_d = GauXC::util::cuda_malloc<double>(rw.size());
  auto* grid_atoms_d = GauXC::util::cuda_malloc<std::int32_t>(grid_atoms.size());
  auto* rho_t_d = GauXC::util::cuda_malloc<double>(rho_t.size());
  auto* gamma_t_d = GauXC::util::cuda_malloc<double>(gamma_t.size());

  auto* omega_d = GauXC::util::cuda_malloc<double>(npts);
  auto* kappa_d = GauXC::util::cuda_malloc<double>(npts);
  auto* domega_drho_d = GauXC::util::cuda_malloc<double>(npts);
  auto* domega_dgamma_d = GauXC::util::cuda_malloc<double>(npts);
  auto* d2omega_drho2_d = GauXC::util::cuda_malloc<double>(npts);
  auto* d2omega_dgamma2_d = GauXC::util::cuda_malloc<double>(npts);
  auto* d2omega_drho_dgamma_d = GauXC::util::cuda_malloc<double>(npts);
  auto* dkappa_drho_d = GauXC::util::cuda_malloc<double>(npts);
  auto* d2kappa_drho2_d = GauXC::util::cuda_malloc<double>(npts);
  auto* F_d = GauXC::util::cuda_malloc<double>(npts);
  auto* U_d = GauXC::util::cuda_malloc<double>(npts);
  auto* W_d = GauXC::util::cuda_malloc<double>(npts);
  auto* eps_d = GauXC::util::cuda_malloc<double>(npts);
  auto* vrho_d = GauXC::util::cuda_malloc<double>(npts);
  auto* vgamma_d = GauXC::util::cuda_malloc<double>(npts);
  auto* grid_gradient_d = GauXC::util::cuda_malloc<double>(3*npts);
  auto* A_d = GauXC::util::cuda_malloc<double>(npts);
  auto* B_d = GauXC::util::cuda_malloc<double>(npts);
  auto* C_d = GauXC::util::cuda_malloc<double>(npts);
  auto* E_d = GauXC::util::cuda_malloc<double>(npts);
  auto* f_rho_t_d = GauXC::util::cuda_malloc<double>(ntrial*npts);
  auto* f_gamma_t_d = GauXC::util::cuda_malloc<double>(ntrial*npts);
  auto* Egr_d = GauXC::util::cuda_malloc<double>(natoms*3*npts);
  auto* Ugr_d = GauXC::util::cuda_malloc<double>(natoms*3*npts);
  auto* Wgr_d = GauXC::util::cuda_malloc<double>(natoms*3*npts);

  GauXC::util::cuda_copy(coords.size(), coords_d, coords.data());
  GauXC::util::cuda_copy(weights.size(), weights_d, weights.data());
  GauXC::util::cuda_copy(rho.size(), rho_d, rho.data());
  GauXC::util::cuda_copy(density_gradient_norm.size(), gamma_d, density_gradient_norm.data());
  GauXC::util::cuda_copy(rw.size(), rw_d, rw.data());
  GauXC::util::cuda_copy(grid_atoms.size(), grid_atoms_d, grid_atoms.data());
  GauXC::util::cuda_copy(rho_t.size(), rho_t_d, rho_t.data());
  GauXC::util::cuda_copy(gamma_t.size(), gamma_t_d, gamma_t.data());

  cudaStream_t stream = 0;
  GauXC::detail::vv10::eval_omega_kappa_cuda(stream, npts, rho_d, gamma_d, params,
    omega_d, kappa_d, domega_drho_d, domega_dgamma_d, d2omega_drho2_d,
    d2omega_dgamma2_d, d2omega_drho_dgamma_d, dkappa_drho_d, d2kappa_drho2_d);
  GauXC::detail::vv10::eval_fock_cuda(stream, npts, coords_d, omega_d, kappa_d,
    npts, coords_d, rw_d, omega_d, kappa_d, F_d, U_d, W_d);
  GauXC::detail::vv10::eval_exc_vxc_from_fock_cuda(stream, npts, rho_d, gamma_d,
    F_d, U_d, W_d, params, eps_d, vrho_d, vgamma_d);
  GauXC::detail::vv10::eval_grid_gradient_cuda(stream, npts, coords_d, omega_d,
    kappa_d, npts, coords_d, rw_d, omega_d, kappa_d, grid_gradient_d);
  GauXC::detail::vv10::eval_hessian_intermediates_cuda(stream, npts, coords_d,
    weights_d, rho_d, omega_d, kappa_d, U_d, W_d, A_d, B_d, C_d, E_d);
  GauXC::detail::vv10::eval_hessian_contraction_cuda(stream, npts, ntrial,
    coords_d, weights_d, rho_d, omega_d, kappa_d, U_d, W_d, A_d, B_d, C_d,
    domega_drho_d, domega_dgamma_d, dkappa_drho_d, d2omega_drho2_d,
    d2omega_dgamma2_d, d2omega_drho_dgamma_d, d2kappa_drho2_d, rho_t_d,
    gamma_t_d, f_rho_t_d, f_gamma_t_d);
  GauXC::detail::vv10::eval_grid_response_cuda(stream, npts, natoms, coords_d,
    weights_d, rho_d, omega_d, kappa_d, grid_atoms_d, Egr_d, Ugr_d, Wgr_d);
  GauXC::util::cuda_device_sync();

  std::vector<double> eps(npts), vrho(npts), vgamma(npts), U(npts), W(npts), E(npts);
  std::vector<double> grid_gradient(3*npts), grid_gradient_ref(3*npts);
  std::vector<double> f_rho_t(ntrial*npts), f_gamma_t(ntrial*npts), Egr(natoms*3*npts);
  GauXC::util::cuda_copy(npts, eps.data(), eps_d);
  GauXC::util::cuda_copy(npts, vrho.data(), vrho_d);
  GauXC::util::cuda_copy(npts, vgamma.data(), vgamma_d);
  GauXC::util::cuda_copy(npts, U.data(), U_d);
  GauXC::util::cuda_copy(npts, W.data(), W_d);
  GauXC::util::cuda_copy(npts, E.data(), E_d);
  GauXC::util::cuda_copy(3*npts, grid_gradient.data(), grid_gradient_d);
  GauXC::util::cuda_copy(ntrial*npts, f_rho_t.data(), f_rho_t_d);
  GauXC::util::cuda_copy(ntrial*npts, f_gamma_t.data(), f_gamma_t_d);
  GauXC::util::cuda_copy(natoms*3*npts, Egr.data(), Egr_d);

  check_vector(eps, {0.00445301813026348, 0.00444920198423763, 0.00445616308451667, 0.00445820144472496, 0.00445560260974276});
  check_vector(vrho, {0.00442962367367506, 0.00442417961353026, 0.00443444131349190, 0.00443937740576438, 0.00443445791849522});
  check_vector(vgamma, {5.4703231981970949e-10, 4.3982859960426448e-09, 3.7195930885603543e-10, 2.0026596148858508e-08, 2.6546741400669848e-09}, 1e-17);
  check_vector(E, {-1.0365967086190920e-04, -1.1129196291360637e-04, -9.7369762355538814e-05, -9.3293041938961706e-05, -9.8490711903360671e-05});
  std::vector<double> omega(npts), kappa(npts);
  GauXC::detail::vv10::eval_omega_kappa(npts, rho.data(), density_gradient_norm.data(), params,
    omega.data(), kappa.data());
  GauXC::detail::vv10::eval_grid_gradient(npts, coords.data(), omega.data(),
    kappa.data(), npts, coords.data(), rw.data(), omega.data(), kappa.data(),
    grid_gradient_ref.data());
  check_vector(grid_gradient, grid_gradient_ref);
  check_vector(f_rho_t, {1.6178037606933446e-07, -1.3249866968720368e-06, 1.1116302821324279e-07, -7.5930158108261938e-08, -1.8068421774630219e-06, -1.0423414119308089e-06, 1.8870904199228002e-06, -1.0109141134429277e-06, 3.5175453056695153e-07, 3.5357534022654595e-07});
  CHECK(Egr[0] == Approx(1.0310153822420043e-05).margin(1e-13));
  CHECK(Egr[29] == Approx(-5.9820261926690516e-06).margin(1e-13));

  GauXC::util::cuda_free(coords_d, weights_d, rho_d, gamma_d, rw_d, grid_atoms_d,
    rho_t_d, gamma_t_d, omega_d, kappa_d, domega_drho_d, domega_dgamma_d,
    d2omega_drho2_d, d2omega_dgamma2_d, d2omega_drho_dgamma_d, dkappa_drho_d,
    d2kappa_drho2_d, F_d, U_d, W_d, eps_d, vrho_d, vgamma_d, grid_gradient_d,
    A_d, B_d, C_d, E_d, f_rho_t_d, f_gamma_t_d, Egr_d, Ugr_d, Wgr_d);
}

TEST_CASE("VV10 CUDA EXC/VXC view helper matches CPU helper", "[vv10][nlc][cuda]") {
  GauXC::IntegratorSettingsNLC settings;
  settings.vv10_b = 6.3;
  settings.vv10_c = 0.0093;
  settings.vv10_tol = 1e-8;

  std::vector<double> eps_cpu(npts, 0.0), vrho_cpu(npts, 0.0), vgamma_cpu(npts, 0.0);
  std::vector<double> eps_cuda(npts, 0.0), vrho_cuda(npts, 0.0), vgamma_cuda(npts, 0.0);

  GauXC::detail::vv10::GridView grid{
    static_cast<std::size_t>(npts), coords.data(), weights.data(), rho.data(), density_gradient_norm.data()
  };

  GauXC::detail::vv10::eval_exc_vxc(
    settings, grid, GauXC::detail::vv10::CorrectionsView{ eps_cpu.data(), vrho_cpu.data(), vgamma_cpu.data() } );
  GauXC::detail::vv10::eval_exc_vxc_cuda(
    0, settings, grid, GauXC::detail::vv10::CorrectionsView{ eps_cuda.data(), vrho_cuda.data(), vgamma_cuda.data() } );

  check_vector(eps_cuda, eps_cpu);
  check_vector(vrho_cuda, vrho_cpu);
  check_vector(vgamma_cuda, vgamma_cpu, 1e-17);
}

TEST_CASE("VV10 CUDA RKS response view helper matches CPU helper", "[vv10][nlc][cuda]") {
  GauXC::IntegratorSettingsNLC settings;
  settings.vv10_b = 6.3;
  settings.vv10_c = 0.0093;
  settings.vv10_tol = 1e-8;

  const std::vector<double> grad_x = {0.08, -0.02, 0.11, 0.04, -0.05};
  const std::vector<double> grad_y = {-0.04, 0.07, 0.01, -0.09, 0.06};
  const std::vector<double> grad_z = {0.03, 0.05, -0.06, 0.02, 0.08};
  const std::vector<double> trial_rho = {0.03, -0.02, 0.04, 0.01, -0.05};
  const std::vector<double> trial_grad_x = {0.002, -0.001, 0.003, -0.002, 0.001};
  const std::vector<double> trial_grad_y = {-0.004, 0.002, 0.001, -0.003, 0.005};
  const std::vector<double> trial_grad_z = {0.001, -0.003, 0.002, 0.004, -0.002};

  std::vector<double> A_cpu(npts, 0.0), B_cpu(3*npts, 0.0);
  std::vector<double> A_cuda(npts, 0.0), B_cuda(3*npts, 0.0);

  GauXC::detail::vv10::ResponseGridView grid{
    static_cast<std::size_t>(npts), coords.data(), weights.data(), rho.data(),
    density_gradient_norm.data(), grad_x.data(), grad_y.data(), grad_z.data()
  };
  GauXC::detail::vv10::TrialView trial{
    trial_rho.data(), trial_grad_x.data(), trial_grad_y.data(), trial_grad_z.data()
  };

  GauXC::detail::vv10::eval_fxc_rks(
    settings, grid, trial, GauXC::detail::vv10::ResponseCorrectionsView{A_cpu.data(), B_cpu.data()} );
  GauXC::detail::vv10::eval_fxc_rks_cuda(
    0, settings, grid, trial, GauXC::detail::vv10::ResponseCorrectionsView{A_cuda.data(), B_cuda.data()} );

  check_vector(A_cuda, A_cpu);
  check_vector(B_cuda, B_cpu);
}

TEST_CASE("VV10 CUDA grid-gradient view helper matches CPU helper", "[vv10][nlc][cuda]") {
  GauXC::IntegratorSettingsNLC settings;
  settings.vv10_b = 6.3;
  settings.vv10_c = 0.0093;
  settings.vv10_tol = 1e-8;

  std::vector<double> gx_cpu(npts, 0.0), gy_cpu(npts, 0.0), gz_cpu(npts, 0.0);
  std::vector<double> gx_cuda(npts, 0.0), gy_cuda(npts, 0.0), gz_cuda(npts, 0.0);

  GauXC::detail::vv10::GridView grid{
    static_cast<std::size_t>(npts), coords.data(), weights.data(), rho.data(), density_gradient_norm.data()
  };

  GauXC::detail::vv10::eval_grid_gradient(
    settings, grid, GauXC::detail::vv10::GridGradientView{gx_cpu.data(), gy_cpu.data(), gz_cpu.data()} );
  GauXC::detail::vv10::eval_grid_gradient_cuda(
    0, settings, grid, GauXC::detail::vv10::GridGradientView{gx_cuda.data(), gy_cuda.data(), gz_cuda.data()} );

  check_vector(gx_cuda, gx_cpu);
  check_vector(gy_cuda, gy_cpu);
  check_vector(gz_cuda, gz_cpu);
}
#endif