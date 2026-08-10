/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy).
 *
 * (c) 2024-2025, Microsoft Corporation
 *
 * All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include <gauxc/grid_factory.hpp>

#include <integratorxx/quadratures/s2/lebedev_laikov.hpp>
#include <integratorxx/quadratures/radial/muraknowles.hpp>
#include <integratorxx/quadratures/radial/mhl.hpp>
#include <integratorxx/quadratures/radial/treutlerahlrichs.hpp>
#include <integratorxx/quadratures/radial/becke.hpp>
#include <integratorxx/composite_quadratures/spherical_quadrature.hpp>
#include <gauxc/exceptions.hpp>

#include <array>
#include <algorithm>

namespace GauXC {

double pyscf_slater_radius_64(AtomicNumber);
double pyscf_gill_radius_93(AtomicNumber);

namespace detail {

using ll_type = IntegratorXX::LebedevLaikov<double>;
using ll_traits = IntegratorXX::quadrature_traits<ll_type>;

template <typename RadialType>
std::vector<double> radial_points_for_spec(UnprunedAtomicGridSpecification unp) {
  RadialType rq(unp.radial_size.get(), unp.radial_scale.get());
  return rq.points();
}

std::vector<double> radial_points_for_spec(UnprunedAtomicGridSpecification unp) {
  switch( unp.radial_quad ) {
    case RadialQuad::Becke:
      return radial_points_for_spec<IntegratorXX::Becke<double,double>>(unp);
    case RadialQuad::MuraKnowles:
      return radial_points_for_spec<IntegratorXX::MuraKnowles<double,double>>(unp);
    case RadialQuad::MurrayHandyLaming:
      return radial_points_for_spec<IntegratorXX::MurrayHandyLaming<double,double>>(unp);
    case RadialQuad::TreutlerAhlrichs:
      return radial_points_for_spec<IntegratorXX::TreutlerAhlrichs<double,double>>(unp);
    default:
      GAUXC_GENERIC_EXCEPTION("Unsupported Radial Quadrature");
      abort();
  }
}

PrunedAtomicGridSpecification make_pruned_spec(
  UnprunedAtomicGridSpecification unp, const std::vector<int64_t>& angular_sizes) {

  if( angular_sizes.size() != static_cast<size_t>(unp.radial_size.get()) )
    GAUXC_GENERIC_EXCEPTION("Invalid PySCF pruning specification");

  std::vector<PruningRegion> pruning_regions;
  size_t idx_st = 0;
  for( size_t i = 1; i <= angular_sizes.size(); ++i ) {
    if( i == angular_sizes.size() || angular_sizes[i] != angular_sizes[idx_st] ) {
      pruning_regions.push_back({idx_st, i, AngularSize(angular_sizes[idx_st])});
      idx_st = i;
    }
  }

  return PrunedAtomicGridSpecification{
    unp.radial_quad, unp.radial_size, unp.radial_scale, pruning_regions
  };
}

const std::array<double, 4>& pyscf_alphas(AtomicNumber Z) {
  static const std::array<double, 4> h_he = {{0.25, 0.5, 1.0, 4.5}};
  static const std::array<double, 4> li_ne = {{0.1667, 0.5, 0.9, 3.5}};
  static const std::array<double, 4> rest = {{0.1, 0.4, 0.8, 2.5}};

  if( Z.get() <= 2 ) return h_he;
  if( Z.get() <= 10 ) return li_ne;
  return rest;
}

std::vector<size_t> pyscf_alpha_regions(
  AtomicNumber Z, UnprunedAtomicGridSpecification unp, double atom_radius ) {

  const auto rads = radial_points_for_spec(unp);
  const auto& alphas = pyscf_alphas(Z);
  std::vector<size_t> regions(rads.size());

  for( size_t i = 0; i < rads.size(); ++i ) {
    size_t region = 0;
    const auto scaled_rad = rads[i] / atom_radius;
    for( const auto alpha : alphas )
      if( scaled_rad > alpha ) ++region;
    regions[i] = region;
  }

  return regions;
}

int64_t previous_lebedev_npts(int64_t npts) {
  const auto order = ll_traits::algebraic_order_by_npts(npts);
  if( order < 0 )
    GAUXC_GENERIC_EXCEPTION("Unsupported PySCF Lebedev angular grid");

  int64_t previous_order = -1;
  auto current_order = ll_traits::next_algebraic_order(3);
  while( current_order < order ) {
    previous_order = current_order;
    current_order = ll_traits::next_algebraic_order(current_order + 1);
  }

  if( previous_order < 0 )
    GAUXC_GENERIC_EXCEPTION("Unsupported PySCF Lebedev angular grid");

  return ll_traits::npts_by_algebraic_order(previous_order);
}

double pyscf_sg1_radius(AtomicNumber Z) {
  const auto radius = pyscf_gill_radius_93(Z);
  if( radius <= 0. )
    GAUXC_GENERIC_EXCEPTION("Atomic number is outside the PySCF SG1 radius table");
  return radius;
}

double pyscf_bragg_radius(AtomicNumber Z) {
  const auto radius = pyscf_slater_radius_64(Z);
  if( radius <= 0. )
    GAUXC_GENERIC_EXCEPTION("Atomic number is outside the PySCF Bragg radius table");
  return radius;
}

} // namespace detail

/*****************/
/**** Visitor ****/
/*****************/
Grid AtomicGridFactory::generate_grid( atomic_grid_variant gs, BatchSize bsz ) {
  return std::visit( [=](auto&& s){ return generate_grid(s, bsz); }, gs );
}

/************************/
/**** Unpruned Grids ****/
/************************/

Grid AtomicGridFactory::generate_unpruned_grid( RadialQuad rq, RadialSize nrad, 
  AngularSize nang, RadialScale rscal, BatchSize bsz) {

  using bk_type  = IntegratorXX::Becke<double, double>;
  using mk_type  = IntegratorXX::MuraKnowles<double,double>;
  using mhl_type = IntegratorXX::MurrayHandyLaming<double,double>;
  using ta_type  = IntegratorXX::TreutlerAhlrichs<double,double>;
  using ll_type  = IntegratorXX::LebedevLaikov<double>;

  ll_type ang_quad( nang.get() );

  switch( rq ) {
    case RadialQuad::Becke:
      return generate_unpruned_grid( bk_type(nrad.get(), rscal.get()),
        std::move(ang_quad), bsz );

    case RadialQuad::MuraKnowles:
      return generate_unpruned_grid( mk_type(nrad.get(), rscal.get()),
        std::move(ang_quad), bsz );

    case RadialQuad::MurrayHandyLaming:
      return generate_unpruned_grid( mhl_type(nrad.get(), rscal.get()),
        std::move(ang_quad), bsz );

    case RadialQuad::TreutlerAhlrichs:
      return generate_unpruned_grid( ta_type(nrad.get(), rscal.get()),
        std::move(ang_quad), bsz );

    default:
      GAUXC_GENERIC_EXCEPTION("Unsupported Radial Quadrature");
      abort();

  }

}

Grid AtomicGridFactory::generate_grid( UnprunedAtomicGridSpecification gs, BatchSize bsz ) {
  return generate_unpruned_grid( gs.radial_quad, gs.radial_size, gs.angular_size,
    gs.radial_scale, bsz );
}




/**********************/
/**** Pruned Grids ****/
/**********************/

template <typename RadialQuad, 
  typename AngularQuad = IntegratorXX::LebedevLaikov<double>>
auto make_pruned_grid(RadialSize nrad, 
  const std::vector<PruningRegion>& pruning_regions,
  RadialScale rscal ) {

  RadialQuad rq(nrad.get(), rscal.get());
  IntegratorXX::RadialGridPartition<AngularQuad> rgp;
  for( auto& region : pruning_regions ) {
    rgp.add_quad( rq, region.idx_st, 
      AngularQuad(region.angular_size.get()) );
  }
  rgp.finalize(rq);

  return std::make_tuple( rq, rgp );

}

Grid AtomicGridFactory::generate_pruned_grid( RadialQuad rq, 
  RadialSize nrad, const std::vector<PruningRegion>& pruning_regions, 
  RadialScale rscal, BatchSize bsz) {

  using mk_type  = IntegratorXX::MuraKnowles<double,double>;
  using mhl_type = IntegratorXX::MurrayHandyLaming<double,double>;
  using ta_type  = IntegratorXX::TreutlerAhlrichs<double,double>;

  switch( rq ) {

    case RadialQuad::MuraKnowles:
    {
      auto [rg, rgp] = 
        make_pruned_grid<mk_type>( nrad, pruning_regions, rscal );
      return generate_pruned_grid(std::move(rg), std::move(rgp), bsz);
    }

    case RadialQuad::MurrayHandyLaming:
    {
      auto [rg, rgp] = 
        make_pruned_grid<mhl_type>( nrad, pruning_regions, rscal );
      return generate_pruned_grid(std::move(rg), std::move(rgp), bsz);
    }

    case RadialQuad::TreutlerAhlrichs:
    {
      auto [rg, rgp] = 
        make_pruned_grid<ta_type>( nrad, pruning_regions, rscal );
      return generate_pruned_grid(std::move(rg), std::move(rgp), bsz);
    }
    case RadialQuad::Becke:
    {
      auto[rg, rgp] = 
        make_pruned_grid<IntegratorXX::Becke<double,double>>( nrad, pruning_regions, rscal );
        return generate_pruned_grid(std::move(rg), std::move(rgp), bsz);
    }

    default:
      GAUXC_GENERIC_EXCEPTION("Unsupported Radial Quadrature");
      abort();

  }

}

Grid AtomicGridFactory::generate_grid( PrunedAtomicGridSpecification gs, BatchSize bsz ) {
  return generate_pruned_grid( gs.radial_quad, gs.radial_size, 
    gs.pruning_regions, gs.radial_scale, bsz );
}



/***************************
 * Default Pruning Schemes *
 ***************************/


PrunedAtomicGridSpecification robust_psi4_pruning_scheme(
  UnprunedAtomicGridSpecification unp ) {

  // Look up order
  // XXX: THIS ONLY WORKS FOR LEBEDEV
  using angular_type = IntegratorXX::LebedevLaikov<double>; 
  using traits = IntegratorXX::quadrature_traits<angular_type>;
  const auto asz = unp.angular_size.get();
  const auto base_order = traits::algebraic_order_by_npts(asz);
  if( base_order < 0 ) GAUXC_GENERIC_EXCEPTION("Invalid Base Grid");

  const auto med_order = 
    traits::next_algebraic_order(base_order > 6 ? base_order-6 : base_order);
  const auto low_order = 7;

  AngularSize med_sz(traits::npts_by_algebraic_order(med_order));
  AngularSize low_sz(traits::npts_by_algebraic_order(low_order));

  // Create Pruning Regions
  const size_t rsz = unp.radial_size.get();
  const size_t r_div_4 = rsz / 4ul + 1ul;
  const size_t r_div_2 = rsz / 2ul + 1ul;
  std::vector<PruningRegion> pruning_regions = {
    {0ul,     r_div_4, low_sz},
    {r_div_4, r_div_2, med_sz},
    {r_div_2,     rsz, unp.angular_size}
  };

  return PrunedAtomicGridSpecification{
    unp.radial_quad, unp.radial_size, unp.radial_scale, pruning_regions
  };
  
}



PrunedAtomicGridSpecification treutler_pruning_scheme(
  UnprunedAtomicGridSpecification unp ) {

  const size_t med_order = 11;
  const size_t low_order = 7;

  // Look up order
  // XXX: THIS ONLY WORKS FOR LEBEDEV
  using angular_type = IntegratorXX::LebedevLaikov<double>;
  using traits = IntegratorXX::quadrature_traits<angular_type>;

  AngularSize med_sz(traits::npts_by_algebraic_order(med_order));
  AngularSize low_sz(traits::npts_by_algebraic_order(low_order));

  // Create Pruning Regions
  const size_t rsz = unp.radial_size.get();
  const size_t r_div_3 = rsz / 3ul + 1ul;
  const size_t r_div_2 = rsz / 2ul + 1ul;
  std::vector<PruningRegion> pruning_regions = {
    {0ul,     r_div_3, low_sz},
    {r_div_3, r_div_2, med_sz},
    {r_div_2, rsz,     unp.angular_size}
  };

  return PrunedAtomicGridSpecification{
    unp.radial_quad, unp.radial_size, unp.radial_scale, pruning_regions
  };
  
}

PrunedAtomicGridSpecification pyscf_sg1_pruning_scheme(
  UnprunedAtomicGridSpecification unp, AtomicNumber Z ) {

  static const std::array<int64_t, 5> leb_ngrid = {{6, 38, 86, 194, 86}};
  const auto regions = detail::pyscf_alpha_regions(
    Z, unp, detail::pyscf_sg1_radius(Z)
  );

  std::vector<int64_t> angular_sizes;
  angular_sizes.reserve(regions.size());
  for( const auto region : regions ) angular_sizes.push_back(leb_ngrid[region]);

  return detail::make_pruned_spec(unp, angular_sizes);
}

PrunedAtomicGridSpecification pyscf_nwchem_pruning_scheme(
  UnprunedAtomicGridSpecification unp, AtomicNumber Z ) {

  const auto n_ang = unp.angular_size.get();
  if( n_ang < 50 ) {
    return detail::make_pruned_spec(
      unp, std::vector<int64_t>(unp.radial_size.get(), n_ang)
    );
  }

  std::array<int64_t, 5> leb_npts;
  if( n_ang == 50 ) {
    leb_npts = {{
      detail::ll_traits::npts_by_algebraic_order(11),
      detail::ll_traits::npts_by_algebraic_order(13),
      detail::ll_traits::npts_by_algebraic_order(13),
      detail::ll_traits::npts_by_algebraic_order(13),
      detail::ll_traits::npts_by_algebraic_order(11)
    }};
  } else {
    const auto previous_n_ang = detail::previous_lebedev_npts(n_ang);
    leb_npts = {{
      detail::ll_traits::npts_by_algebraic_order(11),
      detail::ll_traits::npts_by_algebraic_order(15),
      previous_n_ang,
      n_ang,
      previous_n_ang
    }};
  }

  const auto regions = detail::pyscf_alpha_regions(
    Z, unp, detail::pyscf_bragg_radius(Z)
  );

  std::vector<int64_t> angular_sizes;
  angular_sizes.reserve(regions.size());
  for( const auto region : regions )
    angular_sizes.push_back(leb_npts[region]);

  return detail::make_pruned_spec(unp, angular_sizes);
}

PrunedAtomicGridSpecification pyscf_sgx_pruning_scheme(
  UnprunedAtomicGridSpecification unp, AtomicNumber Z ) {

  static const std::array<std::array<int64_t, 5>, 8> sgx_ang_mapping = {{
    {{ 5,  7, 11, 11,  7}},
    {{ 5,  7, 11, 17, 11}},
    {{ 7, 11, 17, 23, 17}},
    {{ 7, 17, 23, 29, 23}},
    {{ 7, 23, 29, 35, 29}},
    {{11, 29, 35, 41, 35}},
    {{17, 35, 41, 47, 41}},
    {{47, 47, 47, 47, 47}}
  }};

  size_t level = 0;
  const auto n_ang = unp.angular_size.get();
  for( const auto& mapping : sgx_ang_mapping ) {
    const auto level_npts = detail::ll_traits::npts_by_algebraic_order(mapping[3]);
    if( n_ang > level_npts ) ++level;
  }

  if( level >= sgx_ang_mapping.size() )
    GAUXC_GENERIC_EXCEPTION("PySCF SGX pruning is not defined for this angular grid");

  const auto regions = detail::pyscf_alpha_regions(
    Z, unp, detail::pyscf_bragg_radius(Z)
  );

  std::vector<int64_t> angular_sizes;
  angular_sizes.reserve(regions.size());
  for( const auto region : regions )
    angular_sizes.push_back(
      detail::ll_traits::npts_by_algebraic_order(sgx_ang_mapping[level][region])
    );

  return detail::make_pruned_spec(unp, angular_sizes);
}

PrunedAtomicGridSpecification pyscf_treutler_pruning_scheme(
  UnprunedAtomicGridSpecification unp ) {

  using angular_type = IntegratorXX::LebedevLaikov<double>;
  using traits = IntegratorXX::quadrature_traits<angular_type>;

  const size_t rsz = unp.radial_size.get();
  const size_t r_div_3 = rsz / 3ul;
  const size_t r_div_2 = rsz / 2ul;
  std::vector<int64_t> angular_sizes(rsz, unp.angular_size.get());
  std::fill(angular_sizes.begin(), angular_sizes.begin() + r_div_3,
            traits::npts_by_algebraic_order(5));
  std::fill(angular_sizes.begin() + r_div_3, angular_sizes.begin() + r_div_2,
            traits::npts_by_algebraic_order(11));

  return detail::make_pruned_spec(unp, angular_sizes);
}


PrunedAtomicGridSpecification create_pruned_spec(
  PruningScheme scheme, UnprunedAtomicGridSpecification unp
) {

  switch(scheme) {
    case PruningScheme::Robust:
      return robust_psi4_pruning_scheme(unp);
    case PruningScheme::Treutler:
      return treutler_pruning_scheme(unp);
    case PruningScheme::PySCF_Treutler:
      return pyscf_treutler_pruning_scheme(unp);
    case PruningScheme::PySCF_SG1:
      GAUXC_GENERIC_EXCEPTION("PySCF SG1 pruning requires an atomic number");
    case PruningScheme::PySCF_NWChem:
      GAUXC_GENERIC_EXCEPTION("PySCF NWChem pruning requires an atomic number");
    case PruningScheme::PySCF_SGX:
      GAUXC_GENERIC_EXCEPTION("PySCF SGX pruning requires an atomic number");
    
    // Default to Unpruned Grid
    case PruningScheme::Unpruned:
    default:
      std::vector<PruningRegion> pruning_regions = {
        {0ul, (size_t)unp.radial_size.get(), unp.angular_size}
      };
      return PrunedAtomicGridSpecification{
        unp.radial_quad, unp.radial_size, unp.radial_scale, pruning_regions
      };
  }

}

PrunedAtomicGridSpecification create_pruned_spec(
  PruningScheme scheme, AtomicNumber Z, UnprunedAtomicGridSpecification unp
) {

  switch(scheme) {
    case PruningScheme::PySCF_SG1:
      return pyscf_sg1_pruning_scheme(unp, Z);
    case PruningScheme::PySCF_NWChem:
      return pyscf_nwchem_pruning_scheme(unp, Z);
    case PruningScheme::PySCF_SGX:
      return pyscf_sgx_pruning_scheme(unp, Z);
    default:
      return create_pruned_spec(scheme, unp);
  }

}

}
