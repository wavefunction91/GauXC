#pragma once
#include <gauxc/grid.hpp>
#include <integratorxx/composite_quadratures/spherical_quadrature.hpp>
#include <integratorxx/composite_quadratures/pruned_spherical_quadrature.hpp>

namespace GauXC {

struct AtomicGridFactory {

  template <typename RadialType, typename AngularType>
  using unpruned_sphere_type = typename
    IntegratorXX::SphericalQuadrature< std::decay_t<RadialType>, std::decay_t<AngularType>>;
  template <typename RadialType, typename AngularType>
  using pruned_sphere_type = typename
    IntegratorXX::PrunedSphericalQuadrature< std::decay_t<RadialType>, std::decay_t<AngularType>>;

  template <typename RadialType, typename AngularType>
  static Grid generate_unpruned_grid( RadialType&& rq, AngularType&& aq ) {
    using sphere_type = unpruned_sphere_type<RadialType,AngularType>;
    return Grid( std::make_shared<sphere_type>( 
      std::forward<RadialType>(rq), std::forward<AngularType>(aq) 
      )
    );
  }

  static Grid generate_unpruned_grid( RadialQuad, RadialSize, AngularSize, 
                                      RadialScale );

  template <typename RadialType, typename RadialPartitionType>
  static Grid generate_pruned_grid( RadialType&& rq, RadialPartitionType&& rgp ) {
    using angular_type = typename RadialPartitionType::angular_type; 
    using sphere_type = pruned_sphere_type<RadialType,angular_type>;
    return Grid( std::make_shared<sphere_type>( 
      std::forward<RadialType>(rq), std::forward<RadialPartitionType>(rgp)
      )
    );
  }

  template <typename RadialType, typename AngularType>
  static Grid generate_pruned_grid( RadialType&& rq, AngularType&& ang_hgh,
                                    AngularType&& ang_med, AngularType&& ang_low ) {

    const size_t r4 = rq.npts()/4 + !!(rq.npts()%4);
    IntegratorXX::RadialGridPartition 
      rgp( rq, 0,      std::forward<AngularType>(ang_low),
               r4+1,   std::forward<AngularType>(ang_med),
               2*r4+1, std::forward<AngularType>(ang_hgh) );


    return generate_pruned_grid( std::forward<RadialType>(rq), std::move(rgp) );

  }
};

}
