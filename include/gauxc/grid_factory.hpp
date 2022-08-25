#pragma once
#include <gauxc/grid.hpp>
#include <integratorxx/composite_quadratures/spherical_quadrature.hpp>
#include <integratorxx/composite_quadratures/pruned_spherical_quadrature.hpp>

#include <variant>

namespace GauXC {

struct UnprunedAtomicGridSpecification {
  RadialQuad  radial_quad;
  RadialSize  radial_size;
  RadialScale radial_scale;

  AngularSize angular_size;
};

struct PruningRegion {
  size_t idx_st;
  size_t idx_en;
  AngularSize angular_size;

  bool operator==(const PruningRegion& other) const {
    return other.idx_st == idx_st and 
           other.idx_en == idx_en and
           other.angular_size == angular_size;
  }
};

struct PrunedAtomicGridSpecification {
  RadialQuad  radial_quad;
  RadialSize  radial_size;
  RadialScale radial_scale;

  std::vector<PruningRegion> pruning_regions;
};


PrunedAtomicGridSpecification robust_psi4_pruning_scheme(
  UnprunedAtomicGridSpecification
);
PrunedAtomicGridSpecification treutler_pruning_scheme(
  UnprunedAtomicGridSpecification
);

enum class PruningScheme {
  Unpruned,
  Robust,
  Treutler
};

PrunedAtomicGridSpecification create_pruned_spec(
  PruningScheme, UnprunedAtomicGridSpecification
);

using atomic_grid_variant = 
  std::variant<UnprunedAtomicGridSpecification,
               PrunedAtomicGridSpecification>;

struct AtomicGridFactory {

  template <typename RadialType, typename AngularType>
  using unpruned_sphere_type = typename
    IntegratorXX::SphericalQuadrature< 
      std::decay_t<RadialType>, std::decay_t<AngularType>
    >;
  template <typename RadialType, typename AngularType>
  using pruned_sphere_type = typename
    IntegratorXX::PrunedSphericalQuadrature< 
      std::decay_t<RadialType>, std::decay_t<AngularType>
    >;

  template <typename RadialType, typename AngularType>
  static Grid generate_unpruned_grid( RadialType&& rq, AngularType&& aq ) {
    using sphere_type = unpruned_sphere_type<RadialType,AngularType>;
    return Grid( std::make_shared<sphere_type>( 
      std::forward<RadialType>(rq), std::forward<AngularType>(aq) 
      )
    );
  }

  template <typename RadialType, typename RadialPartitionType>
  static Grid generate_pruned_grid( RadialType&& rq, 
    RadialPartitionType&& rgp ) {
    using angular_type = typename std::decay_t<RadialPartitionType>::angular_type; 
    using sphere_type = pruned_sphere_type<RadialType,angular_type>;
    return Grid( std::make_shared<sphere_type>( 
      std::forward<RadialType>(rq), std::forward<RadialPartitionType>(rgp)
      )
    );
  }


  static Grid generate_unpruned_grid( RadialQuad, RadialSize, AngularSize, 
                                      RadialScale );
  static Grid generate_pruned_grid( RadialQuad, RadialSize, 
    const std::vector<PruningRegion>&, RadialScale );


  static Grid generate_grid( UnprunedAtomicGridSpecification gs ); 
  static Grid generate_grid( PrunedAtomicGridSpecification gs ); 

  static Grid generate_grid( atomic_grid_variant gs );

};

}
