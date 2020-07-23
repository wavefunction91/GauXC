#include "catch2/catch.hpp"
#include <gauxc/grid.hpp>

#include <integratorxx/quadratures/lebedev_laikov.hpp>
#include <integratorxx/quadratures/muraknowles.hpp>
#include <integratorxx/quadratures/mhl.hpp>
#include <integratorxx/composite_quadratures/spherical_quadrature.hpp>

#include <random>

using namespace GauXC;

TEST_CASE("Grid", "[grid]") {

  std::default_random_engine gen;
  std::uniform_real_distribution<> pos_real_dist( 0., 10. );

  int64_t n_rad    = 99;
  int64_t n_ang    = 770;
  int64_t batch_sz = 128;

  double r_scal = pos_real_dist(gen);

  IntegratorXX::LebedevLaikov<double>            ang_quad( n_ang         );
  IntegratorXX::MuraKnowles<double,double>       mk_quad ( n_rad, r_scal );
  IntegratorXX::MurrayHandyLaming<double,double> mhl_quad( n_rad, r_scal );

  IntegratorXX::SphericalQuadrature mk_sphere ( mk_quad, ang_quad  );
  IntegratorXX::SphericalQuadrature mhl_sphere( mhl_quad, ang_quad );

  auto mk_batch  = IntegratorXX::make_batcher( batch_sz, mk_sphere  );
  auto mhl_batch = IntegratorXX::make_batcher( batch_sz, mhl_sphere );

  SECTION("Full Construction") {

    RadialQuad rquad = RadialQuad::MuraKnowles;
    SECTION("MK")  { rquad = RadialQuad::MuraKnowles;       }
    SECTION("MHL") { rquad = RadialQuad::MurrayHandyLaming; }

    Grid grid( rquad, RadialSize(n_rad), AngularSize(n_ang), RadialScale(r_scal),
               BatchSize(batch_sz) );


    CHECK( grid.n_rad() == n_rad );
    CHECK( grid.n_ang() == n_ang );
    CHECK( grid.rscal_factor() == r_scal );
    CHECK( grid.max_batch_sz() == batch_sz );

    CHECK( grid.radial_quad() == rquad );

    if( rquad == RadialQuad::MuraKnowles ) {

      for( auto i = 0; i < mk_batch.nbatches(); ++i ) {

        auto&& [box_lo_ref, box_up_ref, points_ref, weights_ref] = mk_batch.at(i);
        auto&& [box_lo, box_up, points, weights] = grid.batcher().at(i);

        CHECK( box_lo_ref == box_lo );
        CHECK( box_up_ref == box_up );

        CHECK( points_ref  == points  );
        CHECK( weights_ref == weights );

      }

    } else if( rquad == RadialQuad::MurrayHandyLaming ) {

      for( auto i = 0; i < mhl_batch.nbatches(); ++i ) {

        auto&& [box_lo_ref, box_up_ref, points_ref, weights_ref] = mhl_batch.at(i);
        auto&& [box_lo, box_up, points, weights] = grid.batcher().at(i);

        CHECK( box_lo_ref == box_lo );
        CHECK( box_up_ref == box_up );

        CHECK( points_ref  == points  );
        CHECK( weights_ref == weights );

      }

    }


    SECTION("Default Batch Size") {
      Grid grid( rquad, RadialSize(n_rad), AngularSize(n_ang), 
                 RadialScale(r_scal) );

      CHECK( grid.max_batch_sz() == 512ll );
    }

    SECTION("Default Radial Quadrature") {
      Grid grid( (RadialSize(n_rad)), AngularSize(n_ang), RadialScale(r_scal),
                 BatchSize(batch_sz) );

      CHECK( grid.radial_quad() == RadialQuad::MuraKnowles );
    }

    SECTION("Default Radial Quadrature + Batch Size") {
      Grid grid( (RadialSize(n_rad)), AngularSize(n_ang), RadialScale(r_scal) ); 
      CHECK( grid.radial_quad() == RadialQuad::MuraKnowles );
      CHECK( grid.max_batch_sz() == 512ll );
    }

  }

}
