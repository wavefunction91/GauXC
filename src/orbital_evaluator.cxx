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
#include <gauxc/orbital_evaluator.hpp>

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <gauxc/exceptions.hpp>
#include <gauxc/xc_integrator/local_work_driver.hpp>

#include "xc_integrator/local_work_driver/host/blas.hpp"
#include "xc_integrator/local_work_driver/host/local_host_work_driver.hpp"

namespace GauXC {

namespace {

/// Per-thread point batch size, sized so the AO scratch buffer
/// (nbf*batch doubles) fits comfortably in L2 (~16 MB target).
inline int64_t choose_batch_size(int32_t nbf) {
  constexpr int64_t kTargetBytesPerThread = 16ll * 1024 * 1024;
  constexpr int64_t kMinBatch = 256;
  constexpr int64_t kMaxBatch = 8192;
  if (nbf <= 0) return kMaxBatch;
  const int64_t b =
      kTargetBytesPerThread / (static_cast<int64_t>(sizeof(double)) * nbf);
  return std::clamp(b, kMinBatch, kMaxBatch);
}

/// Single-block submat_map for the no-screening case (nbe == nbf); makes
/// `eval_xmat` route directly to GEMM without invoking the gather path.
inline LocalHostWorkDriver::submat_map_t full_submat_map(int32_t nbf) {
  return {{ {int32_t{0}, nbf, int32_t{0}} }};
}

}  // namespace

struct OrbitalEvaluator::Impl {
  BasisSet<double> basis;
  std::unique_ptr<LocalWorkDriver> driver_owner;
  LocalHostWorkDriver* host_driver = nullptr;  // non-owning view
  std::vector<int32_t> shell_list;
  int32_t nbf_ = 0;

  void init(BasisSet<double> bs, ExecutionSpace exec) {
    if (exec != ExecutionSpace::Host) {
      GAUXC_GENERIC_EXCEPTION(
          "OrbitalEvaluator: only ExecutionSpace::Host is currently supported.");
    }

    basis = std::move(bs);

    driver_owner = LocalWorkDriverFactory::make_local_work_driver(
        ExecutionSpace::Host, "Reference");
    host_driver = dynamic_cast<LocalHostWorkDriver*>(driver_owner.get());
    if (host_driver == nullptr) {
      GAUXC_GENERIC_EXCEPTION(
          "OrbitalEvaluator: LocalWorkDriverFactory did not return a "
          "LocalHostWorkDriver.");
    }

    shell_list.resize(basis.size());
    std::iota(shell_list.begin(), shell_list.end(), int32_t{0});
    nbf_ = basis.nbf();
  }
};

OrbitalEvaluator::OrbitalEvaluator(BasisSet<double> basis, ExecutionSpace exec)
    : pimpl_(std::make_unique<Impl>()) {
  pimpl_->init(std::move(basis), exec);
}

OrbitalEvaluator::~OrbitalEvaluator() noexcept = default;
OrbitalEvaluator::OrbitalEvaluator(OrbitalEvaluator&&) noexcept = default;
OrbitalEvaluator& OrbitalEvaluator::operator=(OrbitalEvaluator&&) noexcept =
    default;

int32_t OrbitalEvaluator::nbf() const noexcept { return pimpl_->nbf_; }

const BasisSet<double>& OrbitalEvaluator::basis() const noexcept {
  return pimpl_->basis;
}

void OrbitalEvaluator::eval_orbital(int64_t npts, const double* points,
                                    const double* C, double* out) const {
  eval_orbitals(npts, points, /*nmo=*/1, C, /*ldc=*/pimpl_->nbf_, out,
                /*ldo=*/npts);
}

void OrbitalEvaluator::eval_orbitals(int64_t npts, const double* points,
                                     int32_t nmo, const double* C, int64_t ldc,
                                     double* out, int64_t ldo) const {
  if (npts == 0 || nmo == 0) return;
  if (points == nullptr || C == nullptr || out == nullptr) {
    GAUXC_GENERIC_EXCEPTION(
        "OrbitalEvaluator::eval_orbitals: null pointer argument.");
  }
  const int32_t nbf = pimpl_->nbf_;
  if (ldc < nbf) {
    GAUXC_GENERIC_EXCEPTION(
        "OrbitalEvaluator::eval_orbitals: ldc must be >= nbf().");
  }
  if (ldo < npts) {
    GAUXC_GENERIC_EXCEPTION(
        "OrbitalEvaluator::eval_orbitals: ldo must be >= npts.");
  }

  const int32_t nshells = pimpl_->basis.nshells();
  const int64_t batch_size = choose_batch_size(nbf);
  const int64_t n_batches = (npts + batch_size - 1) / batch_size;
  LocalHostWorkDriver* driver = pimpl_->host_driver;
  const int32_t* shell_list = pimpl_->shell_list.data();
  const BasisSet<double>& basis = pimpl_->basis;

#pragma omp parallel
  {
    std::vector<double> ao_buf(static_cast<size_t>(nbf) * batch_size);

#pragma omp for schedule(dynamic, 1)
    for (int64_t b = 0; b < n_batches; ++b) {
      const int64_t p0 = b * batch_size;
      const int64_t np = std::min<int64_t>(batch_size, npts - p0);

      driver->eval_collocation(static_cast<size_t>(np),
                               static_cast<size_t>(nshells),
                               static_cast<size_t>(nbf), points + 3 * p0,
                               basis, shell_list, ao_buf.data());

      // out_slab(np, nmo) = ao^T(np, nbf) @ C(nbf, nmo); j-th MO column
      // lives at out + j*ldo + p0.
      blas::gemm<double>(
          /*TA=*/'T', /*TB=*/'N',
          /*M=*/static_cast<int>(np), /*N=*/static_cast<int>(nmo),
          /*K=*/nbf, /*ALPHA=*/1.0, /*A=*/ao_buf.data(), /*LDA=*/nbf,
          /*B=*/C, /*LDB=*/static_cast<int>(ldc), /*BETA=*/0.0,
          /*C=*/out + p0, /*LDC=*/static_cast<int>(ldo));
    }
  }
}

void OrbitalEvaluator::eval_density(int64_t npts, const double* points,
                                    const double* D, int64_t ldd,
                                    double* out) const {
  if (npts == 0) return;
  if (points == nullptr || D == nullptr || out == nullptr) {
    GAUXC_GENERIC_EXCEPTION(
        "OrbitalEvaluator::eval_density: null pointer argument.");
  }
  const int32_t nbf = pimpl_->nbf_;
  if (ldd < nbf) {
    GAUXC_GENERIC_EXCEPTION(
        "OrbitalEvaluator::eval_density: ldd must be >= nbf().");
  }

  const int32_t nshells = pimpl_->basis.nshells();
  const int64_t batch_size = choose_batch_size(nbf);
  const int64_t n_batches = (npts + batch_size - 1) / batch_size;
  LocalHostWorkDriver* driver = pimpl_->host_driver;
  const int32_t* shell_list = pimpl_->shell_list.data();
  const BasisSet<double>& basis = pimpl_->basis;

#pragma omp parallel
  {
    std::vector<double> ao_buf(static_cast<size_t>(nbf) * batch_size);
    std::vector<double> dm_ao_buf(static_cast<size_t>(nbf) * batch_size);
    // Sized for the eval_xmat submat-gather contract; unused on the
    // nbe == nbf fast path but allocated unconditionally for safety.
    std::vector<double> xmat_scr(static_cast<size_t>(nbf) * nbf);
    auto submat_map = full_submat_map(nbf);

#pragma omp for schedule(dynamic, 1)
    for (int64_t b = 0; b < n_batches; ++b) {
      const int64_t p0 = b * batch_size;
      const int64_t np = std::min<int64_t>(batch_size, npts - p0);

      driver->eval_collocation(static_cast<size_t>(np),
                               static_cast<size_t>(nshells),
                               static_cast<size_t>(nbf), points + 3 * p0,
                               basis, shell_list, ao_buf.data());

      // dm_ao = D @ ao
      driver->eval_xmat(
          /*npts=*/static_cast<size_t>(np), /*nbf=*/static_cast<size_t>(nbf),
          /*nbe=*/static_cast<size_t>(nbf), submat_map, /*fac=*/1.0,
          /*P=*/D, /*ldp=*/static_cast<size_t>(ldd),
          /*basis_eval=*/ao_buf.data(), /*ldb=*/static_cast<size_t>(nbf),
          /*X=*/dm_ao_buf.data(), /*ldx=*/static_cast<size_t>(nbf),
          /*scr=*/xmat_scr.data());

      // rho[p] = sum_mu ao(mu, p) * dm_ao(mu, p)
      driver->eval_uvvar_lda_rks(
          /*npts=*/static_cast<size_t>(np), /*nbe=*/static_cast<size_t>(nbf),
          /*basis_eval=*/ao_buf.data(),
          /*X=*/dm_ao_buf.data(), /*ldx=*/static_cast<size_t>(nbf),
          /*den_eval=*/out + p0);
    }
  }
}

}  // namespace GauXC
