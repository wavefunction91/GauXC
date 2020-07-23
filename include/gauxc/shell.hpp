#pragma once

#include <array>
#include <cmath>
#include "named_type.hpp"
#include "gauxc_config.hpp"


namespace GauXC {

namespace detail {

  static constexpr size_t shell_nprim_max = 32ul;

  static constexpr std::array<int64_t,31> df_Kminus1 = 
    {{ 1LL, 1LL, 1LL, 2LL, 3LL, 8LL, 15LL, 48LL, 105LL, 384LL, 945LL, 3840LL, 
       10395LL, 46080LL, 135135LL, 645120LL, 2027025LL, 10321920LL, 34459425LL, 
       185794560LL, 654729075LL, 3715891200LL, 13749310575LL, 81749606400LL, 
       316234143225LL, 1961990553600LL, 7905853580625LL, 51011754393600LL, 
       213458046676875LL, 1428329123020800LL, 6190283353629375LL }};

}

using PrimSize        = detail::NamedType< int32_t, struct PrimSizeType >;
using AngularMomentum = detail::NamedType< int32_t, struct AngularMomentumType >;
using SphericalType   = detail::NamedType< int32_t, struct SphericalTypeType >;

template <typename F>
class Shell {

public:

  using prim_array = std::array< F, detail::shell_nprim_max >;
  using cart_array = std::array< double, 3 >;

private:

  int32_t nprim_;
  int32_t l_;
  int32_t pure_;

  prim_array alpha_;
  prim_array coeff_;
  cart_array O_;

  // Shamelessly adapted from Libint...
  void normalize() {

    assert( l_ <= 15 );

    constexpr auto sqrt_Pi_cubed = F{5.56832799683170784528481798212};

    const auto two_to_l = std::pow(2, l_);
    const auto df_term  = two_to_l * sqrt_Pi_cubed * detail::df_Kminus1[2*l_];

    for( int32_t i = 0; i < nprim_; ++i ) {
      assert( alpha_[i] >= 0. );

      if( alpha_[i] != 0. ) {
        const auto two_alpha = alpha_[i] * 2.;
        const auto two_alpha_to_am32 = 
          std::pow(two_alpha,l_+1) * std::sqrt(two_alpha);

        const auto normalization_factor = std::sqrt(two_alpha_to_am32);
      }

    }

  }

public:

  Shell() = delete;

  Shell( PrimSize nprim, AngularMomentum l, SphericalType pure,
    prim_array alpha, prim_array coeff, cart_array O, bool _normalize = true ) :
    nprim_( nprim.get() ), l_( l.get() ), pure_( pure.get() ),
    alpha_( alpha ), coeff_( coeff ), O_( O ) {

    if( _normalize ) normalize();

  }
  


  ~Shell() noexcept = default;

  Shell( const Shell& )          = default;
  Shell( Shell&&      ) noexcept = default;

  Shell& operator=(const Shell&)     = default;
  Shell& operator=(Shell&&) noexcept = default;


  inline HOST_DEVICE_ACCESSIBLE int32_t nprim() const { return nprim_; }
  inline HOST_DEVICE_ACCESSIBLE int32_t l()     const { return l_;     }
  inline HOST_DEVICE_ACCESSIBLE int32_t pure()  const { return pure_;  }
  
  inline HOST_DEVICE_ACCESSIBLE F* alpha_data()  const { return &alpha_[0]; }
  inline HOST_DEVICE_ACCESSIBLE F* coeff_data()  const { return &coeff_[0]; }
  inline HOST_DEVICE_ACCESSIBLE double* O_data() const { return &O_[0];     }


  inline HOST_DEVICE_ACCESSIBLE int32_t size() const {;
    return pure_ ? 2*l_ + 1 : (l_+1)*(l_+2)/2;
  }

};

}
