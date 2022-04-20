#pragma once
#include <gauxc/shell.hpp>
#include <gauxc/basisset.hpp>
#include <gauxc/exceptions.hpp>

namespace GauXC {
namespace detail {
  struct cartesian_point {
    double x, y, z;
  };

  static constexpr size_t nprim_pair_max = 64ul;

  template <typename Integral>
  inline constexpr Integral packed_lt_index(Integral i, Integral j, Integral m) {
    return i + ((2*m - j - 1)*j)/2;
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

  //std::vector< PrimitivePair<F> > prim_pairs_;
  std::array< PrimitivePair<F>, detail::nprim_pair_max > prim_pairs_;
  size_t nprim_pairs_ = 0;

  void generate( const_shell_ref bra, const_shell_ref ket ) {
    //prim_pairs_.resize( bra.nprim() * ket.nprim() );

    detail::cartesian_point A{ bra.O()[0], bra.O()[1], bra.O()[2] };
    detail::cartesian_point B{ ket.O()[0], ket.O()[1], ket.O()[2] };

    const auto rABx = A.x - B.x;
    const auto rABy = A.y - B.y;
    const auto rABz = A.z - B.z;

    const auto dAB = rABx*rABx + rABy*rABy + rABz*rABz;

    const auto np_bra = bra.nprim();
    const auto np_ket = ket.nprim();
    nprim_pairs_ = 0;
    for( auto i = 0, ij = 0; i < np_bra; ++i       )
    for( auto j = 0;         j < np_ket; ++j, ++ij ) {
      if( nprim_pairs_ >= detail::nprim_pair_max ) 
        GAUXC_GENERIC_EXCEPTION("Too Many Primitive Pairs");

      nprim_pairs_++;
      auto& pair = prim_pairs_[ij];
      const auto alpha_bra = bra.alpha()[i];
      const auto alpha_ket = ket.alpha()[j];

      const auto g    = alpha_bra + alpha_ket;
      const auto oo_g = 1 / g;

      pair.P.x = (alpha_bra * A.x + alpha_ket * B.x) * oo_g;
      pair.P.y = (alpha_bra * A.y + alpha_ket * B.y) * oo_g;
      pair.P.z = (alpha_bra * A.z + alpha_ket * B.z) * oo_g;

      pair.PA.x = pair.P.x - A.x;
      pair.PA.y = pair.P.y - A.y;
      pair.PA.z = pair.P.z - A.z;

      pair.PB.x = pair.P.x - B.x;
      pair.PB.y = pair.P.y - B.y;
      pair.PB.z = pair.P.z - B.z;

      pair.K_coeff_prod = 2 * M_PI * oo_g *
        bra.coeff()[i] * ket.coeff()[j] *
        std::exp( -alpha_bra * alpha_ket * dAB * oo_g );
      pair.gamma = g;
      pair.gamma_inv = oo_g;
    } // loop over prim pairs
  } // generate

public:

  ShellPair( const Shell<F>& bra, const Shell<F>& ket ) {
    if( bra.l() >= ket.l() ) generate(bra,ket);
    else                     generate(ket,bra);
  }

  PrimitivePair<F>* prim_pairs() { return prim_pairs_.data(); }
  const PrimitivePair<F>* prim_pairs() const { return prim_pairs_.data(); }

  size_t nprim_pairs() const { return nprim_pairs_; }

};


template <typename F>
class ShellPairCollection {
  size_t nshells_ = 0;
  std::vector<ShellPair<F>> shell_pairs_;

public:
  ShellPairCollection( const BasisSet<F>& basis ) {
    nshells_ = basis.size();

    // Pack into column-major lower triangle
    shell_pairs_.reserve( (nshells_ * (nshells_+1))/2 );
    for(int j = 0; j < nshells_; ++j)
    for(int i = j; i < nshells_; ++i) {
      shell_pairs_.emplace_back( basis[i], basis[j] );
    }
  }

  // Retreive unique LT element
  inline auto& at( size_t i, size_t j ) {
    if( j <= i ) {
      const auto idx = detail::packed_lt_index(i,j,nshells_);
      return shell_pairs_[idx];
    } else return at(j,i);
  }

  inline const auto& at( size_t i, size_t j ) const {
    if( j <= i ) {
      const auto idx = detail::packed_lt_index(i,j,nshells_);
      return shell_pairs_[idx];
    } else return at(j,i);
  }

};

}
