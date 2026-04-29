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
#include <array>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <vector>

#include <omp.h>

#include <gauxc/basisset_map.hpp>
#include <gauxc/external/cube.hpp>
#include <gauxc/exceptions.hpp>
#include <gauxc/molecule.hpp>
#include <gauxc/xc_integrator/local_work_driver.hpp>

#include "xc_integrator/integrator_util/integrator_common.hpp"
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

/// Fill `pts[3*np]` with grid-point coordinates for batch [p0, p0+np).
/// Avoids materialising the full 3*N point array for grid-based evaluation.
inline void fill_batch_pts(const CubeGrid& g, int64_t p0, int64_t np,
                           double* pts) {
  const int64_t nz = g.nz, ny = g.ny;
  for (int64_t i = 0; i < np; ++i) {
    const int64_t k  = p0 + i;
    const int64_t iz = k % nz;
    const int64_t iy = (k / nz) % ny;
    const int64_t ix = k / (ny * nz);
    pts[3 * i + 0] = g.origin[0] + static_cast<double>(ix) * g.spacing[0];
    pts[3 * i + 1] = g.origin[1] + static_cast<double>(iy) * g.spacing[1];
    pts[3 * i + 2] = g.origin[2] + static_cast<double>(iz) * g.spacing[2];
  }
}

/// Axis-aligned bbox of a set of npts AoS points (length 3*npts).
struct PointBbox {
  std::array<double, 3> lo;
  std::array<double, 3> hi;
};
inline PointBbox compute_bbox(const double* points, int64_t npts) {
  PointBbox b{{points[0], points[1], points[2]},
              {points[0], points[1], points[2]}};
  for (int64_t p = 1; p < npts; ++p) {
    const double* xyz = points + 3 * p;
    for (int k = 0; k < 3; ++k) {
      if (xyz[k] < b.lo[k]) b.lo[k] = xyz[k];
      if (xyz[k] > b.hi[k]) b.hi[k] = xyz[k];
    }
  }
  return b;
}

/// Squared distance from `center` to the nearest point in the axis-aligned
/// bbox `[lo, hi]^3`. Zero if the center lies inside the bbox.
inline double dist2_center_to_bbox(const double* center,
                                   const PointBbox& bbox) {
  double d2 = 0.0;
  for (int k = 0; k < 3; ++k) {
    const double c = center[k];
    if (c < bbox.lo[k]) {
      const double dx = bbox.lo[k] - c;
      d2 += dx * dx;
    } else if (c > bbox.hi[k]) {
      const double dx = c - bbox.hi[k];
      d2 += dx * dx;
    }
  }
  return d2;
}

}  // namespace

struct OrbitalEvaluator::Impl {
  BasisSet<double> basis;
  std::unique_ptr<LocalWorkDriver> driver_owner;
  LocalHostWorkDriver* host_driver = nullptr;  // non-owning view
  std::vector<int32_t> shell_list;
  // BasisSetMap powers shell -> AO range lookups for the screened
  // submat_map. We don't have a Molecule available so we pass an empty
  // one; the only field that needs it (shell_to_center) is unused here.
  std::unique_ptr<BasisSetMap> basis_map;
  // Per-shell squared cutoff radius, cached so the screening loop is just
  // a comparison instead of a square root per shell per batch.
  std::vector<double> shell_cutoff_r2;
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

    basis_map = std::make_unique<BasisSetMap>(basis, Molecule{});

    shell_cutoff_r2.resize(basis.size());
    for (size_t s = 0; s < basis.size(); ++s) {
      const double r = basis[s].cutoff_radius();
      shell_cutoff_r2[s] = r * r;
    }
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

  const int32_t nshells_total = pimpl_->basis.nshells();
  const int64_t batch_size = choose_batch_size(nbf);
  const int64_t n_batches = (npts + batch_size - 1) / batch_size;
  LocalHostWorkDriver* driver = pimpl_->host_driver;
  const BasisSet<double>& basis = pimpl_->basis;
  const auto& shell_cutoff_r2 = pimpl_->shell_cutoff_r2;
  const BasisSetMap& basis_map = *pimpl_->basis_map;

#pragma omp parallel
  {
    std::vector<double> ao_buf(static_cast<size_t>(nbf) * batch_size);
    std::vector<double> C_compressed(static_cast<size_t>(nbf) * nmo);
    std::vector<int32_t> screened_shells;
    screened_shells.reserve(nshells_total);

#pragma omp for schedule(dynamic, 1)
    for (int64_t b = 0; b < n_batches; ++b) {
      const int64_t p0 = b * batch_size;
      const int64_t np = std::min<int64_t>(batch_size, npts - p0);

      // Per-batch shell screening: keep shells whose cutoff_radius reaches
      // any point of this batch's bounding box.
      const PointBbox bbox = compute_bbox(points + 3 * p0, np);
      screened_shells.clear();
      int32_t nbe = 0;
      for (int32_t s = 0; s < nshells_total; ++s) {
        if (dist2_center_to_bbox(basis[s].O_data(), bbox) <
            shell_cutoff_r2[s]) {
          screened_shells.push_back(s);
          nbe += basis[s].size();
        }
      }
      if (nbe == 0) {
        // All shells screen out for this batch -> orbital is identically
        // zero. Zero the slab and skip eval.
        for (int32_t j = 0; j < nmo; ++j) {
          double* out_col = out + static_cast<size_t>(j) * ldo + p0;
          std::fill(out_col, out_col + np, 0.0);
        }
        continue;
      }

      driver->eval_collocation(
          static_cast<size_t>(np),
          static_cast<size_t>(screened_shells.size()),
          static_cast<size_t>(nbe), points + 3 * p0, basis,
          screened_shells.data(), ao_buf.data());

      // Gather the rows of C corresponding to surviving shells into a
      // contiguous (nbe, nmo) col-major buffer so the contraction is a
      // dense GEMM.
      {
        int32_t row = 0;
        for (int32_t s : screened_shells) {
          const auto rng = basis_map.shell_to_ao_range(s);
          const int32_t shell_nbe = rng.second - rng.first;
          for (int32_t j = 0; j < nmo; ++j) {
            const double* C_col = C + static_cast<size_t>(j) * ldc;
            double* dst = C_compressed.data() +
                          static_cast<size_t>(j) * nbe + row;
            std::copy(C_col + rng.first, C_col + rng.second, dst);
          }
          row += shell_nbe;
        }
      }

      // out_slab(np, nmo) = ao^T(np, nbe) @ C_compressed(nbe, nmo); j-th MO
      // column lives at out + j*ldo + p0.
      blas::gemm<double>(
          /*TA=*/'T', /*TB=*/'N',
          /*M=*/static_cast<int>(np), /*N=*/static_cast<int>(nmo),
          /*K=*/nbe, /*ALPHA=*/1.0, /*A=*/ao_buf.data(), /*LDA=*/nbe,
          /*B=*/C_compressed.data(), /*LDB=*/nbe, /*BETA=*/0.0,
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

  const int32_t nshells_total = pimpl_->basis.nshells();
  const int64_t batch_size = choose_batch_size(nbf);
  const int64_t n_batches = (npts + batch_size - 1) / batch_size;
  LocalHostWorkDriver* driver = pimpl_->host_driver;
  const BasisSet<double>& basis = pimpl_->basis;
  const auto& shell_cutoff_r2 = pimpl_->shell_cutoff_r2;
  const BasisSetMap& basis_map = *pimpl_->basis_map;

#pragma omp parallel
  {
    std::vector<double> ao_buf(static_cast<size_t>(nbf) * batch_size);
    std::vector<double> dm_ao_buf(static_cast<size_t>(nbf) * batch_size);
    // eval_xmat scratch when the submat-gather path is taken; sized for
    // worst case (nbe == nbf).
    std::vector<double> xmat_scr(static_cast<size_t>(nbf) * nbf);
    std::vector<int32_t> screened_shells;
    screened_shells.reserve(nshells_total);

#pragma omp for schedule(dynamic, 1)
    for (int64_t b = 0; b < n_batches; ++b) {
      const int64_t p0 = b * batch_size;
      const int64_t np = std::min<int64_t>(batch_size, npts - p0);

      const PointBbox bbox = compute_bbox(points + 3 * p0, np);
      screened_shells.clear();
      int32_t nbe = 0;
      for (int32_t s = 0; s < nshells_total; ++s) {
        if (dist2_center_to_bbox(basis[s].O_data(), bbox) <
            shell_cutoff_r2[s]) {
          screened_shells.push_back(s);
          nbe += basis[s].size();
        }
      }
      if (nbe == 0) {
        std::fill(out + p0, out + p0 + np, 0.0);
        continue;
      }

      driver->eval_collocation(
          static_cast<size_t>(np),
          static_cast<size_t>(screened_shells.size()),
          static_cast<size_t>(nbe), points + 3 * p0, basis,
          screened_shells.data(), ao_buf.data());

      // Build the compressed submat_map for this batch (gather the rows
      // of D corresponding to surviving shells). Falls through to a single
      // direct-gemm submat_map when nbe == nbf.
      LocalHostWorkDriver::submat_map_t submat_map;
      if (nbe == nbf) {
        submat_map = full_submat_map(nbf);
      } else {
        std::tie(submat_map, std::ignore) =
            gen_compressed_submat_map(basis_map, screened_shells, nbf, nbf);
      }

      // dm_ao = D_compressed @ ao
      driver->eval_xmat(
          /*npts=*/static_cast<size_t>(np), /*nbf=*/static_cast<size_t>(nbf),
          /*nbe=*/static_cast<size_t>(nbe), submat_map, /*fac=*/1.0,
          /*P=*/D, /*ldp=*/static_cast<size_t>(ldd),
          /*basis_eval=*/ao_buf.data(), /*ldb=*/static_cast<size_t>(nbe),
          /*X=*/dm_ao_buf.data(), /*ldx=*/static_cast<size_t>(nbe),
          /*scr=*/xmat_scr.data());

      // rho[p] = sum_mu ao(mu, p) * dm_ao(mu, p)
      driver->eval_uvvar_lda_rks(
          /*npts=*/static_cast<size_t>(np), /*nbe=*/static_cast<size_t>(nbe),
          /*basis_eval=*/ao_buf.data(),
          /*X=*/dm_ao_buf.data(), /*ldx=*/static_cast<size_t>(nbe),
          /*den_eval=*/out + p0);
    }
  }
}

// ---------------------------------------------------------------------------
// CubeGrid overloads (3b): generate per-batch point coordinates on-the-fly,
// avoiding the 3*num_points()*8 byte temporary coordinate array.
// ---------------------------------------------------------------------------

void OrbitalEvaluator::eval_orbital(const CubeGrid& grid,
                                    const double* C, double* out) const {
  eval_orbitals(grid, /*nmo=*/1, C, /*ldc=*/pimpl_->nbf_, out,
                /*ldo=*/grid.num_points());
}

void OrbitalEvaluator::eval_orbitals(const CubeGrid& grid,
                                     int32_t nmo, const double* C, int64_t ldc,
                                     double* out, int64_t ldo) const {
  const int64_t npts = grid.num_points();
  if (npts == 0 || nmo == 0) return;
  if (C == nullptr || out == nullptr) {
    GAUXC_GENERIC_EXCEPTION(
        "OrbitalEvaluator::eval_orbitals(grid): null pointer argument.");
  }
  const int32_t nbf = pimpl_->nbf_;
  if (ldc < nbf) {
    GAUXC_GENERIC_EXCEPTION(
        "OrbitalEvaluator::eval_orbitals(grid): ldc must be >= nbf().");
  }
  if (ldo < npts) {
    GAUXC_GENERIC_EXCEPTION(
        "OrbitalEvaluator::eval_orbitals(grid): ldo must be >= npts.");
  }

  const int32_t nshells_total = pimpl_->basis.nshells();
  const int64_t batch_size = choose_batch_size(nbf);
  const int64_t n_batches = (npts + batch_size - 1) / batch_size;
  LocalHostWorkDriver* driver = pimpl_->host_driver;
  const BasisSet<double>& basis = pimpl_->basis;
  const auto& shell_cutoff_r2 = pimpl_->shell_cutoff_r2;
  const BasisSetMap& basis_map = *pimpl_->basis_map;

  // Choose between task pipeline and omp-for based on saturation.
  // When n_batches >> nthreads, omp-for already fully saturates all cores
  // and the task pipeline's scheduling overhead is a net negative.
  const int nthreads = omp_get_max_threads();
  const bool use_pipeline = (n_batches < 2 * nthreads);

  if (use_pipeline) {
  // Pipelined evaluation (3d): overlap collocation of batch B+1 with
  // GEMM of batch B using OMP tasks with dependency tags.
  struct Slot {
    std::vector<double> pts;
    std::vector<double> ao;
    std::vector<double> C_comp;
    std::vector<int32_t> shells;
    int32_t nbe;
    int64_t p0, np;
  };

  const int nslots =
      std::min<int64_t>(nthreads, n_batches);
  std::vector<Slot> slots(nslots);
  for (auto& s : slots) {
    s.pts.resize(static_cast<size_t>(batch_size) * 3);
    s.ao.resize(static_cast<size_t>(nbf) * batch_size);
    s.C_comp.resize(static_cast<size_t>(nbf) * nmo);
    s.shells.reserve(nshells_total);
    s.nbe = 0;
    s.p0 = 0;
    s.np = 0;
  }

  // Dependency tag array — raw array for omp depend().
  std::unique_ptr<char[]> dtag(new char[nslots]{});
  char* dtag_ptr = dtag.get();
  (void)dtag_ptr;  // used only in omp depend() clauses

#pragma omp parallel
#pragma omp single
  {
    for (int64_t b = 0; b < n_batches; ++b) {
      const int si = static_cast<int>(b % nslots);

      // Phase A: generate points, screen shells, eval collocation.
#pragma omp task shared(slots, grid, basis, shell_cutoff_r2) \
    depend(out: dtag_ptr[si]) firstprivate(b, si)
      {
        Slot& sl = slots[si];
        sl.p0 = b * batch_size;
        sl.np = std::min<int64_t>(batch_size, npts - sl.p0);

        fill_batch_pts(grid, sl.p0, sl.np, sl.pts.data());
        const PointBbox bbox = compute_bbox(sl.pts.data(), sl.np);
        sl.shells.clear();
        sl.nbe = 0;
        for (int32_t s = 0; s < nshells_total; ++s) {
          if (dist2_center_to_bbox(basis[s].O_data(), bbox) <
              shell_cutoff_r2[s]) {
            sl.shells.push_back(s);
            sl.nbe += basis[s].size();
          }
        }
        if (sl.nbe > 0) {
          driver->eval_collocation(
              static_cast<size_t>(sl.np),
              static_cast<size_t>(sl.shells.size()),
              static_cast<size_t>(sl.nbe), sl.pts.data(), basis,
              sl.shells.data(), sl.ao.data());
        }
      }

      // Phase B: gather C + GEMM → out.
#pragma omp task shared(slots, out, C, basis_map) \
    depend(in: dtag_ptr[si]) firstprivate(b, si)
      {
        Slot& sl = slots[si];
        if (sl.nbe == 0) {
          for (int32_t j = 0; j < nmo; ++j) {
            double* oc = out + static_cast<size_t>(j) * ldo + sl.p0;
            std::fill(oc, oc + sl.np, 0.0);
          }
        } else {
          const int32_t nbe = sl.nbe;
          {
            int32_t row = 0;
            for (int32_t s : sl.shells) {
              const auto rng = basis_map.shell_to_ao_range(s);
              const int32_t shell_nbe = rng.second - rng.first;
              for (int32_t j = 0; j < nmo; ++j) {
                const double* C_col = C + static_cast<size_t>(j) * ldc;
                double* dst = sl.C_comp.data() +
                              static_cast<size_t>(j) * nbe + row;
                std::copy(C_col + rng.first, C_col + rng.second, dst);
              }
              row += shell_nbe;
            }
          }

          blas::gemm<double>(
              'T', 'N', static_cast<int>(sl.np), static_cast<int>(nmo),
              nbe, 1.0, sl.ao.data(), nbe, sl.C_comp.data(), nbe,
              0.0, out + sl.p0, static_cast<int>(ldo));
        }
      }
    }
  }

  } else {
    // Saturated path: plain omp-for, lower overhead at high thread counts.

#pragma omp parallel
  {
    std::vector<double> batch_pts(static_cast<size_t>(batch_size) * 3);
    std::vector<double> ao_buf(static_cast<size_t>(nbf) * batch_size);
    std::vector<double> C_compressed(static_cast<size_t>(nbf) * nmo);
    std::vector<int32_t> screened_shells;
    screened_shells.reserve(nshells_total);

#pragma omp for schedule(dynamic, 1)
    for (int64_t b = 0; b < n_batches; ++b) {
      const int64_t p0 = b * batch_size;
      const int64_t np = std::min<int64_t>(batch_size, npts - p0);

      fill_batch_pts(grid, p0, np, batch_pts.data());
      const PointBbox bbox = compute_bbox(batch_pts.data(), np);
      screened_shells.clear();
      int32_t nbe = 0;
      for (int32_t s = 0; s < nshells_total; ++s) {
        if (dist2_center_to_bbox(basis[s].O_data(), bbox) <
            shell_cutoff_r2[s]) {
          screened_shells.push_back(s);
          nbe += basis[s].size();
        }
      }
      if (nbe == 0) {
        for (int32_t j = 0; j < nmo; ++j) {
          double* out_col = out + static_cast<size_t>(j) * ldo + p0;
          std::fill(out_col, out_col + np, 0.0);
        }
        continue;
      }

      driver->eval_collocation(
          static_cast<size_t>(np),
          static_cast<size_t>(screened_shells.size()),
          static_cast<size_t>(nbe), batch_pts.data(), basis,
          screened_shells.data(), ao_buf.data());

      {
        int32_t row = 0;
        for (int32_t s : screened_shells) {
          const auto rng = basis_map.shell_to_ao_range(s);
          const int32_t shell_nbe = rng.second - rng.first;
          for (int32_t j = 0; j < nmo; ++j) {
            const double* C_col = C + static_cast<size_t>(j) * ldc;
            double* dst = C_compressed.data() +
                          static_cast<size_t>(j) * nbe + row;
            std::copy(C_col + rng.first, C_col + rng.second, dst);
          }
          row += shell_nbe;
        }
      }

      blas::gemm<double>(
          'T', 'N', static_cast<int>(np), static_cast<int>(nmo), nbe,
          1.0, ao_buf.data(), nbe, C_compressed.data(), nbe,
          0.0, out + p0, static_cast<int>(ldo));
    }
  }
  }  // end saturated fallback
}

void OrbitalEvaluator::eval_density(const CubeGrid& grid,
                                    const double* D, int64_t ldd,
                                    double* out) const {
  const int64_t npts = grid.num_points();
  if (npts == 0) return;
  if (D == nullptr || out == nullptr) {
    GAUXC_GENERIC_EXCEPTION(
        "OrbitalEvaluator::eval_density(grid): null pointer argument.");
  }
  const int32_t nbf = pimpl_->nbf_;
  if (ldd < nbf) {
    GAUXC_GENERIC_EXCEPTION(
        "OrbitalEvaluator::eval_density(grid): ldd must be >= nbf().");
  }

  const int32_t nshells_total = pimpl_->basis.nshells();
  const int64_t batch_size = choose_batch_size(nbf);
  const int64_t n_batches = (npts + batch_size - 1) / batch_size;
  LocalHostWorkDriver* driver = pimpl_->host_driver;
  const BasisSet<double>& basis = pimpl_->basis;
  const auto& shell_cutoff_r2 = pimpl_->shell_cutoff_r2;
  const BasisSetMap& basis_map = *pimpl_->basis_map;

#pragma omp parallel
  {
    std::vector<double> batch_pts(static_cast<size_t>(batch_size) * 3);
    std::vector<double> ao_buf(static_cast<size_t>(nbf) * batch_size);
    std::vector<double> dm_ao_buf(static_cast<size_t>(nbf) * batch_size);
    std::vector<double> xmat_scr(static_cast<size_t>(nbf) * nbf);
    std::vector<int32_t> screened_shells;
    screened_shells.reserve(nshells_total);

#pragma omp for schedule(dynamic, 1)
    for (int64_t b = 0; b < n_batches; ++b) {
      const int64_t p0 = b * batch_size;
      const int64_t np = std::min<int64_t>(batch_size, npts - p0);

      fill_batch_pts(grid, p0, np, batch_pts.data());
      const PointBbox bbox = compute_bbox(batch_pts.data(), np);
      screened_shells.clear();
      int32_t nbe = 0;
      for (int32_t s = 0; s < nshells_total; ++s) {
        if (dist2_center_to_bbox(basis[s].O_data(), bbox) <
            shell_cutoff_r2[s]) {
          screened_shells.push_back(s);
          nbe += basis[s].size();
        }
      }
      if (nbe == 0) {
        std::fill(out + p0, out + p0 + np, 0.0);
        continue;
      }

      driver->eval_collocation(
          static_cast<size_t>(np),
          static_cast<size_t>(screened_shells.size()),
          static_cast<size_t>(nbe), batch_pts.data(), basis,
          screened_shells.data(), ao_buf.data());

      LocalHostWorkDriver::submat_map_t submat_map;
      if (nbe == nbf) {
        submat_map = full_submat_map(nbf);
      } else {
        std::tie(submat_map, std::ignore) =
            gen_compressed_submat_map(basis_map, screened_shells, nbf, nbf);
      }

      driver->eval_xmat(
          static_cast<size_t>(np), static_cast<size_t>(nbf),
          static_cast<size_t>(nbe), submat_map, 1.0,
          D, static_cast<size_t>(ldd),
          ao_buf.data(), static_cast<size_t>(nbe),
          dm_ao_buf.data(), static_cast<size_t>(nbe),
          xmat_scr.data());

      driver->eval_uvvar_lda_rks(
          static_cast<size_t>(np), static_cast<size_t>(nbe),
          ao_buf.data(),
          dm_ao_buf.data(), static_cast<size_t>(nbe),
          out + p0);
    }
  }
}

}  // namespace GauXC
