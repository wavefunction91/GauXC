#include "integral_bounds.hpp"
#include "rys_integral.h"
#include <vector>
#include <gauxc/util/geometry.hpp>


namespace GauXC {
namespace util  {

struct RysShell {
  shells shell;
  RysShell( int nprim, int l ) {
    shell.m = nprim;
    shell.L = l;

    shell.coeff = new coefficients[nprim];
  }

  RysShell( const Shell<double>& sh ) :
    RysShell( sh.nprim(), sh.l() ) {

    for( auto i = 0; i < shell.m; ++i ) {
      shell.coeff[i].coeff = sh.coeff()[i];
      shell.coeff[i].alpha = sh.alpha()[i];
    }

  }

  ~RysShell() noexcept {
    if( shell.coeff ) delete shell.coeff;
  }
};


template <int n, typename T>
inline constexpr T integral_pow( T x ) {
  if constexpr (n == 1) return x;
  else                  return x * integral_pow<n-1>(x);
}


inline constexpr double max_coulomb_00( double Kab, double gamma ) {
  return Kab / gamma;
}

inline constexpr double max_coulomb_20( double Kab, double Rab, double alpha, double beta, 
  double gamma ) {
  return Kab * ( gamma + Rab * integral_pow<2>(beta) ) / integral_pow<3>(gamma);
}

inline constexpr double max_coulomb_22( double Kab, double Rab, double alpha, double beta, 
  double gamma ) {
  return Kab / integral_pow<5>(gamma) * 
    ( Rab * integral_pow<3>(alpha) +
      alpha * beta * (4. - Rab*beta ) +
      integral_pow<2>(beta) * (2.  + Rab * beta ) +
      integral_pow<2>(alpha) * (2. - Rab * beta + Rab*Rab * beta*beta)
    );
}

inline double max_coulomb( int l_a, int l_b, double Kab, double Rab, double alpha, 
  double beta, double gamma ) {

  constexpr double pi2 = 2. * M_PI;
  if( l_a == 0 and l_b == 0 ) return pi2 * max_coulomb_00( Kab, gamma );
  if( l_a == 2 and l_b == 2 ) return pi2 * max_coulomb_22( Kab, Rab, alpha, beta, gamma );
  if( l_a == 2 and l_b == 0 ) return pi2 * max_coulomb_20( Kab, Rab, alpha, beta, gamma );
  if( l_a == 0 and l_b == 2 ) return pi2 * max_coulomb_20( Kab, Rab, beta, alpha, gamma );

  const int l_a_p = l_a + (l_a % 2);
  const int l_b_p = l_b + (l_b % 2);

  const int l_a_m = l_a - (l_a % 2);
  const int l_b_m = l_b - (l_b % 2);

  if( l_a_p > 2 or l_b_p > 2 ) throw std::runtime_error("Case Not Handled"); 

  double V_pm = std::numeric_limits<double>::infinity();
  if( l_a_p == 0 and l_b_m == 0 ) 
    V_pm = max_coulomb_00( Kab, gamma );
  else if( l_a_p == 2 and l_b_m == 0 ) 
    V_pm = max_coulomb_20( Kab, Rab, alpha, beta, gamma );
  else if( l_a_p == 0 and l_b_m == 2 ) 
    V_pm = max_coulomb_20( Kab, Rab, beta, alpha, gamma );
  else if( l_a_p == 2 and l_b_m == 2 )
    V_pm = max_coulomb_22( Kab, Rab, alpha, beta, gamma );

  double V_mp = std::numeric_limits<double>::infinity();
  if( l_a_m == 0 and l_b_p == 0 ) 
    V_mp = max_coulomb_00( Kab, gamma );
  else if( l_a_p == 2 and l_b_m == 0 ) 
    V_mp = max_coulomb_20( Kab, Rab, alpha, beta, gamma );
  else if( l_a_m == 0 and l_b_p == 2 ) 
    V_mp = max_coulomb_20( Kab, Rab, beta, alpha, gamma );
  else if( l_a_m == 2 and l_b_p == 2 )
    V_mp = max_coulomb_22( Kab, Rab, alpha, beta, gamma );

  return pi2 * std::sqrt(V_pm * V_mp);
}


template <typename T>
T max_coulomb( const Shell<T>& bra, const Shell<T>& ket) {

  const auto A = bra.O();
  const auto B = ket.O();
  const auto RAB = geometry::euclidean_dist( A, B );

  double max_val = 0.;
  for( auto i = 0; i < bra.nprim(); ++i )
  for( auto j = 0; j < ket.nprim(); ++j ) {
    const auto alpha = bra.alpha()[i];
    const auto beta  = ket.alpha()[j];
    const auto gamma = alpha + beta;

    const auto Kab = std::exp( - alpha*beta*RAB / gamma );

    const auto c_a = bra.coeff()[i];
    const auto c_b = ket.coeff()[j];
    const auto c = std::abs( c_a * c_b );

    max_val += c * max_coulomb( bra.l(), ket.l(), Kab, RAB, alpha, beta, gamma );
  }

  return max_val;
}


template double max_coulomb( const Shell<double>&, const Shell<double>& );

}
}
