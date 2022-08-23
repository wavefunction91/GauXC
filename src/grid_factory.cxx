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

}
