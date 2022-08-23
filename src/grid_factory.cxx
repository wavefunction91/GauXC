#include <gauxc/grid_factory.hpp>

#include <integratorxx/quadratures/lebedev_laikov.hpp>
#include <integratorxx/quadratures/muraknowles.hpp>
#include <integratorxx/quadratures/mhl.hpp>
#include <integratorxx/quadratures/treutleraldrichs.hpp>
#include <integratorxx/composite_quadratures/spherical_quadrature.hpp>
#include <gauxc/exceptions.hpp>

namespace GauXC {

/*****************/
/**** Visitor ****/
/*****************/
Grid AtomicGridFactory::generate_grid( atomic_grid_variant gs ) {
  return std::visit( [](auto&& s){ return generate_grid(s); }, gs );
}

/************************/
/**** Unpruned Grids ****/
/************************/

Grid AtomicGridFactory::generate_unpruned_grid( RadialQuad rq, RadialSize nrad, 
  AngularSize nang, RadialScale rscal) {

  using mk_type  = IntegratorXX::MuraKnowles<double,double>;
  using mhl_type = IntegratorXX::MurrayHandyLaming<double,double>;
  using ta_type  = IntegratorXX::TreutlerAldrichs<double,double>;
  using ll_type  = IntegratorXX::LebedevLaikov<double>;

  ll_type ang_quad( nang.get() );

  switch( rq ) {

    case RadialQuad::MuraKnowles:
      return generate_unpruned_grid( mk_type(nrad.get(), rscal.get()),
        std::move(ang_quad) );

    case RadialQuad::MurrayHandyLaming:
      return generate_unpruned_grid( mhl_type(nrad.get(), rscal.get()),
        std::move(ang_quad) );

    case RadialQuad::TreutlerAldrichs:
      return generate_unpruned_grid( ta_type(nrad.get(), rscal.get()),
        std::move(ang_quad) );

    default:
      GAUXC_GENERIC_EXCEPTION("Unsupported Radial Quadrature");
      abort();

  }

}

Grid AtomicGridFactory::generate_grid( UnprunedAtomicGridSpecification gs ) {
  return generate_unpruned_grid( gs.radial_quad, gs.radial_size, gs.angular_size,
    gs.radial_scale );
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
  RadialScale rscal) {

  using mk_type  = IntegratorXX::MuraKnowles<double,double>;
  using mhl_type = IntegratorXX::MurrayHandyLaming<double,double>;
  using ta_type  = IntegratorXX::TreutlerAldrichs<double,double>;

  switch( rq ) {

    case RadialQuad::MuraKnowles:
    {
      auto [rg, rgp] = 
        make_pruned_grid<mk_type>( nrad, pruning_regions, rscal );
      return generate_pruned_grid(std::move(rg), std::move(rgp));
    }

    case RadialQuad::MurrayHandyLaming:
    {
      auto [rg, rgp] = 
        make_pruned_grid<mhl_type>( nrad, pruning_regions, rscal );
      return generate_pruned_grid(std::move(rg), std::move(rgp));
    }

    case RadialQuad::TreutlerAldrichs:
    {
      auto [rg, rgp] = 
        make_pruned_grid<ta_type>( nrad, pruning_regions, rscal );
      return generate_pruned_grid(std::move(rg), std::move(rgp));
    }

    default:
      GAUXC_GENERIC_EXCEPTION("Unsupported Radial Quadrature");
      abort();

  }

}

Grid AtomicGridFactory::generate_grid( PrunedAtomicGridSpecification gs ) {
  return generate_pruned_grid( gs.radial_quad, gs.radial_size, 
    gs.pruning_regions, gs.radial_scale );
}



/***************************
 * Default Pruning Schemes *
 ***************************/


PrunedAtomicGridSpecification robust_psi4_pruning_scheme(
  UnprunedAtomicGridSpecification unp ) {

  // Look up order
  // XXX: THIS ONLY WORKS FOR LEBEDEV
  using namespace IntegratorXX::detail::lebedev;
  const auto asz = unp.angular_size.get();
  const auto base_order = algebraic_order_by_npts(asz);
  if( base_order < 0 ) GAUXC_GENERIC_EXCEPTION("Invalid Base Grid");

  const auto med_order = 
    next_algebraic_order(base_order > 6 ? base_order-6 : base_order);
  const auto low_order = 7;

  AngularSize med_sz(npts_by_algebraic_order(med_order));
  AngularSize low_sz(npts_by_algebraic_order(low_order));

  // Create Pruning Regions
  const size_t rsz = unp.radial_size.get();
  const size_t r_div_4 = rsz / 4ul;
  std::vector<PruningRegion> pruning_regions = {
    {  0ul,       r_div_4,   low_sz},
    {  r_div_4, 2ul*r_div_4, med_sz},
    {2ul*r_div_4,       rsz, unp.angular_size}
  };

  return PrunedAtomicGridSpecification{
    unp.radial_quad, unp.radial_size, unp.radial_scale, pruning_regions
  };
  
}

}
