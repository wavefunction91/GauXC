/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include <gauxc/shell.hpp>
#include <gauxc/basisset.hpp>
#include <gauxc/exceptions.hpp>

namespace GauXC {
namespace detail {
  struct cartesian_point {
    double x, y, z;
  };

#ifdef GAUXC_SHELLPAIR_PRIMPAIR_MAX
  static constexpr size_t nprim_pair_max = GAUXC_SHELLPAIR_PRIMPAIR_MAX;
#else
  static constexpr size_t nprim_pair_max = 64ul;
#endif

  template <typename Integral>
  inline intmax_t csr_index( size_t i, size_t j, Integral* row_ptr, Integral* col_ind ) {
    const auto j_st = col_ind + row_ptr[i];
    const auto j_en = col_ind + row_ptr[i+1];
    auto it = std::lower_bound(j_st, j_en, j);
    if( it != j_en and *it == j )
      return std::distance(col_ind, it);
    else return -1;

  }
}

template <typename F>
struct PrimitivePair {
  detail::cartesian_point P;
  detail::cartesian_point PA;
  detail::cartesian_point PB;

  F K_coeff_prod;
  F gamma;
  F gamma_inv;
};

template <typename F>
class ShellPair {

  using shell_type = Shell<F>;
  using const_shell_ref = const shell_type&;

  std::array< PrimitivePair<F>, detail::nprim_pair_max > prim_pairs_;
  size_t nprim_pairs_ = 0;

  void generate( const_shell_ref bra, const_shell_ref ket ) {

    detail::cartesian_point A{ bra.O()[0], bra.O()[1], bra.O()[2] };
    detail::cartesian_point B{ ket.O()[0], ket.O()[1], ket.O()[2] };

    const auto rABx = A.x - B.x;
    const auto rABy = A.y - B.y;
    const auto rABz = A.z - B.z;

    const auto dAB = rABx*rABx + rABy*rABy + rABz*rABz;

    const auto np_bra = bra.nprim();
    const auto np_ket = ket.nprim();
    nprim_pairs_ = 0;
    for( auto i = 0; i < np_bra; ++i )
    for( auto j = 0; j < np_ket; ++j ) {
      if( nprim_pairs_ >= detail::nprim_pair_max ) 
        GAUXC_GENERIC_EXCEPTION("Too Many Primitive Pairs");

      auto& pair = prim_pairs_[nprim_pairs_];
      const auto alpha_bra = bra.alpha()[i];
      const auto alpha_ket = ket.alpha()[j];

      const auto g    = alpha_bra + alpha_ket;
      const auto oo_g = 1 / g;

      const auto Kab = 2 * M_PI * oo_g *
        bra.coeff()[i] * ket.coeff()[j] *
        std::exp( -alpha_bra * alpha_ket * dAB * oo_g );

      if(std::abs(Kab) < 1e-12) continue;

      nprim_pairs_++;
      pair.P.x = (alpha_bra * A.x + alpha_ket * B.x) * oo_g;
      pair.P.y = (alpha_bra * A.y + alpha_ket * B.y) * oo_g;
      pair.P.z = (alpha_bra * A.z + alpha_ket * B.z) * oo_g;

      pair.PA.x = pair.P.x - A.x;
      pair.PA.y = pair.P.y - A.y;
      pair.PA.z = pair.P.z - A.z;

      pair.PB.x = pair.P.x - B.x;
      pair.PB.y = pair.P.y - B.y;
      pair.PB.z = pair.P.z - B.z;

      pair.K_coeff_prod = Kab;
      pair.gamma = g;
      pair.gamma_inv = oo_g;
    } // loop over prim pairs
  } // generate

public:

  ShellPair( ) : nprim_pairs_(0) {}

  ShellPair( const Shell<F>& bra, const Shell<F>& ket ) {
    if( bra.l() >= ket.l() ) generate(bra,ket);
    else                     generate(ket,bra);
  }

  inline HOST_DEVICE_ACCESSIBLE PrimitivePair<F>* prim_pairs() { return detail::contiguous_data(prim_pairs_); }
  inline HOST_DEVICE_ACCESSIBLE const PrimitivePair<F>* prim_pairs() const { return detail::contiguous_data(prim_pairs_); }

  inline HOST_DEVICE_ACCESSIBLE size_t nprim_pairs() const { return nprim_pairs_; }

};


template <typename F>
class ShellPairCollection {
  size_t nshells_ = 0;
  std::vector<ShellPair<F>> shell_pairs_;
  std::vector<size_t> row_ptr_, col_ind_;
  ShellPair<F> dummy;

public:
  ShellPairCollection( const BasisSet<F>& basis ) {
    nshells_ = basis.size();

    // Sparse Storage based on primitive screening
    row_ptr_.resize(nshells_+1);
    row_ptr_[0] = 0;
    for(size_t i = 0; i < nshells_; ++i) {

      size_t nnz_row = 0;
      for(size_t j = 0; j <= i; ++j) {
        ShellPair<F> sp(basis[i], basis[j]);
        if(sp.nprim_pairs()) {
          nnz_row++;
          col_ind_.emplace_back(j);
          shell_pairs_.emplace_back(std::move(sp));
        }
      }
      row_ptr_[i+1] = row_ptr_[i] + nnz_row;
    }    
  }

  inline int64_t get_linear_shell_pair_index(size_t i, size_t j) const {
    return detail::csr_index(i, j, row_ptr_.data(), col_ind_.data());
  }

  // Retreive unique LT element
  inline auto& at( size_t i, size_t j ) {
    auto idx = get_linear_shell_pair_index(i,j);
    return idx >= 0 ? shell_pairs_[idx] : dummy;
  }

  inline const auto& at( size_t i, size_t j ) const {
    auto idx = get_linear_shell_pair_index(i,j);
    return idx >= 0 ? shell_pairs_[idx] : dummy;
  }

  inline size_t nshells() const { return nshells_; }
  inline size_t npairs() const { return shell_pairs_.size(); }
  inline auto* shell_pairs() { return shell_pairs_.data(); }
  inline auto* shell_pairs() const { return shell_pairs_.data(); }

  inline auto& row_ptr() { return row_ptr_; }
  inline auto& row_ptr() const { return row_ptr_; }
  inline auto& col_ind() { return col_ind_; }
  inline auto& col_ind() const { return col_ind_; }

};

}
