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
#include <gauxc/external/cube.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <gauxc/exceptions.hpp>

namespace GauXC {

// =============================================================================
// CubeGrid
// =============================================================================

CubeGrid CubeGrid::from_molecule(const Molecule& mol, int64_t nx, int64_t ny,
                                 int64_t nz, double margin) {
  if (mol.empty()) {
    GAUXC_GENERIC_EXCEPTION(
        "CubeGrid::from_molecule: molecule has no atoms.");
  }
  if (nx < 1 || ny < 1 || nz < 1) {
    GAUXC_GENERIC_EXCEPTION(
        "CubeGrid::from_molecule: nx, ny, nz must be >= 1.");
  }

  double xmin = mol[0].x, xmax = mol[0].x;
  double ymin = mol[0].y, ymax = mol[0].y;
  double zmin = mol[0].z, zmax = mol[0].z;
  for (const auto& a : mol) {
    xmin = std::min(xmin, a.x);
    xmax = std::max(xmax, a.x);
    ymin = std::min(ymin, a.y);
    ymax = std::max(ymax, a.y);
    zmin = std::min(zmin, a.z);
    zmax = std::max(zmax, a.z);
  }

  CubeGrid grid;
  grid.origin = {xmin - margin, ymin - margin, zmin - margin};
  grid.nx = nx;
  grid.ny = ny;
  grid.nz = nz;
  const double ex = (xmax - xmin) + 2.0 * margin;
  const double ey = (ymax - ymin) + 2.0 * margin;
  const double ez = (zmax - zmin) + 2.0 * margin;
  grid.spacing[0] = nx > 1 ? ex / static_cast<double>(nx - 1) : 0.0;
  grid.spacing[1] = ny > 1 ? ey / static_cast<double>(ny - 1) : 0.0;
  grid.spacing[2] = nz > 1 ? ez / static_cast<double>(nz - 1) : 0.0;
  return grid;
}

std::vector<double> CubeGrid::points() const {
  std::vector<double> pts(static_cast<size_t>(num_points()) * 3);
  points_into(pts.data());
  return pts;
}

void CubeGrid::points_into(double* out) const {
  size_t k = 0;
  for (int64_t ix = 0; ix < nx; ++ix) {
    const double x = origin[0] + spacing[0] * static_cast<double>(ix);
    for (int64_t iy = 0; iy < ny; ++iy) {
      const double y = origin[1] + spacing[1] * static_cast<double>(iy);
      for (int64_t iz = 0; iz < nz; ++iz) {
        const double z = origin[2] + spacing[2] * static_cast<double>(iz);
        out[3 * k + 0] = x;
        out[3 * k + 1] = y;
        out[3 * k + 2] = z;
        ++k;
      }
    }
  }
}

// Hand-rolled %13.5E formatter. Required field is fixed 13 chars
// (sign + D + '.' + 5d + 'E' + sign + 2d). snprintf in the inner loop
// dominates write time on large grids; this matches glibc snprintf
// bit-for-bit on the fast path (1-2 digit exponents) and falls back to
// snprintf for 3-digit-exponent edge cases. Output buffer must be
// exactly 13 bytes (no trailing NUL).

namespace {

inline void format_e13_5(double v, char* out) {
  if (std::isnan(v)) {
    std::memcpy(out, "          NaN", 13);
    return;
  }
  if (std::isinf(v)) {
    std::memcpy(out, v < 0 ? "         -Inf" : "          Inf", 13);
    return;
  }

  bool negative = std::signbit(v);
  double absv = std::fabs(v);

  if (absv == 0.0) {
    // glibc "%13.5E":  0.0 -> "  0.00000E+00";  -0.0 -> " -0.00000E+00".
    out[0] = ' ';
    out[1] = negative ? '-' : ' ';
    out[2] = '0';
    out[3] = '.';
    out[4] = '0';
    out[5] = '0';
    out[6] = '0';
    out[7] = '0';
    out[8] = '0';
    out[9] = 'E';
    out[10] = '+';
    out[11] = '0';
    out[12] = '0';
    return;
  }

  // Exponent via floor(log10), with corrections for FP edge cases
  // (e.g. 9.99999 rounding up across a power-of-ten boundary).
  int exp10 = static_cast<int>(std::floor(std::log10(absv)));
  double scale = std::pow(10.0, -exp10);
  double mant = absv * scale;

  long long mant_int = static_cast<long long>(std::llround(mant * 1e5));
  if (mant_int >= 1000000) {
    mant_int = 100000;
    ++exp10;
  } else if (mant_int < 100000) {
    --exp10;
    scale = std::pow(10.0, -exp10);
    mant = absv * scale;
    mant_int = static_cast<long long>(std::llround(mant * 1e5));
    if (mant_int >= 1000000) {
      mant_int = 999999;
    } else if (mant_int < 100000) {
      mant_int = 100000;
    }
  }

  // 3+ digit exponents have a different field layout; defer to snprintf.
  if (exp10 > 99 || exp10 < -99) {
    char tmp[32];
    const int n = std::snprintf(tmp, sizeof(tmp), "%13.5E", v);
    if (n >= 13) {
      std::memcpy(out, tmp + (n - 13), 13);
    } else {
      const int pad = 13 - n;
      for (int i = 0; i < pad; ++i) out[i] = ' ';
      std::memcpy(out + pad, tmp, static_cast<size_t>(n));
    }
    return;
  }

  // Fast path. Layout: [sp][sign][D][.][d4 d3 d2 d1 d0][E][esign][e1 e0]
  out[0] = ' ';
  out[1] = negative ? '-' : ' ';

  char digits[6];
  for (int i = 5; i >= 0; --i) {
    digits[i] = static_cast<char>('0' + (mant_int % 10));
    mant_int /= 10;
  }
  out[2] = digits[0];
  out[3] = '.';
  out[4] = digits[1];
  out[5] = digits[2];
  out[6] = digits[3];
  out[7] = digits[4];
  out[8] = digits[5];
  out[9] = 'E';
  out[10] = exp10 < 0 ? '-' : '+';

  const int aexp = exp10 < 0 ? -exp10 : exp10;
  out[11] = static_cast<char>('0' + (aexp / 10));
  out[12] = static_cast<char>('0' + (aexp % 10));
}

}  // namespace

// =============================================================================
// write_cube
// =============================================================================

void write_cube(const std::string& path, const Molecule& mol,
                const CubeGrid& grid, const double* field,
                const std::string& comment) {
  if (field == nullptr) {
    GAUXC_GENERIC_EXCEPTION("write_cube: field pointer is null.");
  }
  if (grid.num_points() <= 0) {
    GAUXC_GENERIC_EXCEPTION("write_cube: grid has zero points.");
  }

  std::FILE* f = std::fopen(path.c_str(), "w");
  if (f == nullptr) {
    GAUXC_GENERIC_EXCEPTION("write_cube: failed to open output file: " + path);
  }

  // --- Header ---
  std::fprintf(f, "%s\n",
               comment.empty() ? "GauXC cube file" : comment.c_str());
  std::fprintf(f, "Generated by GauXC\n");

  // natoms + origin (Bohr).
  std::fprintf(f, "%5lld %12.6f %12.6f %12.6f\n",
               static_cast<long long>(mol.size()), grid.origin[0],
               grid.origin[1], grid.origin[2]);

  // Three voxel-axis lines (axis-aligned grid).
  std::fprintf(f, "%5lld %12.6f %12.6f %12.6f\n",
               static_cast<long long>(grid.nx), grid.spacing[0], 0.0, 0.0);
  std::fprintf(f, "%5lld %12.6f %12.6f %12.6f\n",
               static_cast<long long>(grid.ny), 0.0, grid.spacing[1], 0.0);
  std::fprintf(f, "%5lld %12.6f %12.6f %12.6f\n",
               static_cast<long long>(grid.nz), 0.0, 0.0, grid.spacing[2]);

  // One line per atom: Z, partial charge (0.0), x, y, z (Bohr).
  for (const auto& atom : mol) {
    std::fprintf(f, "%5lld %12.6f %12.6f %12.6f %12.6f\n",
                 static_cast<long long>(atom.Z.get()), 0.0, atom.x, atom.y,
                 atom.z);
  }

  // --- Data block ---
  // Each (ix, iy) row is grid.nz values, six per line, %13.5E. The cube
  // format requires a newline at the end of every (ix, iy) row regardless
  // of how many values land on the last line. Rows are independent, so we
  // format them in parallel into a pre-sized buffer and commit with a
  // single fwrite.
  const int64_t nz = grid.nz;
  const int64_t lines_per_row = (nz + 5) / 6;
  // Worst case 6*13 + 1 = 79 bytes per line; over-allocates the trailing
  // line of each row but avoids a precise sizing pass.
  const int64_t bytes_per_row = lines_per_row * (6 * 13 + 1);
  const int64_t n_rows = grid.nx * grid.ny;
  std::vector<char> buf(static_cast<size_t>(bytes_per_row * n_rows));
  std::vector<int64_t> row_byte_count(static_cast<size_t>(n_rows), 0);

#pragma omp parallel for schedule(static)
  for (int64_t row = 0; row < n_rows; ++row) {
    const double* row_data =
        field + static_cast<size_t>(row) * static_cast<size_t>(nz);
    char* dst = buf.data() + static_cast<size_t>(row) *
                                 static_cast<size_t>(bytes_per_row);
    int64_t off = 0;
    for (int64_t iz = 0; iz < nz; ++iz) {
      format_e13_5(row_data[iz], dst + off);
      off += 13;
      // Every 6 values OR at the end of the row → newline.
      if (((iz + 1) % 6 == 0) || (iz + 1 == nz)) {
        dst[off++] = '\n';
      }
    }
    row_byte_count[static_cast<size_t>(row)] = off;
  }

  // Compact rows in place (worst-case padding between them) and emit
  // with a single fwrite.
  if (n_rows > 1) {
    int64_t write_off = row_byte_count[0];
    for (int64_t row = 1; row < n_rows; ++row) {
      const int64_t src_off = row * bytes_per_row;
      const int64_t len = row_byte_count[static_cast<size_t>(row)];
      std::memmove(buf.data() + write_off, buf.data() + src_off,
                   static_cast<size_t>(len));
      write_off += len;
    }
    if (std::fwrite(buf.data(), 1, static_cast<size_t>(write_off), f) !=
        static_cast<size_t>(write_off)) {
      std::fclose(f);
      GAUXC_GENERIC_EXCEPTION("write_cube: short write to " + path);
    }
  } else {
    if (std::fwrite(buf.data(), 1, static_cast<size_t>(row_byte_count[0]),
                    f) != static_cast<size_t>(row_byte_count[0])) {
      std::fclose(f);
      GAUXC_GENERIC_EXCEPTION("write_cube: short write to " + path);
    }
  }

  if (std::fclose(f) != 0) {
    GAUXC_GENERIC_EXCEPTION("write_cube: failed to close " + path);
  }
}

}  // namespace GauXC
