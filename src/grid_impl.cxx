#include "grid_impl.hpp"

#include <integratorxx/quadratures/lebedev_laikov.hpp>
#include <integratorxx/quadratures/muraknowles.hpp>
#include <integratorxx/quadratures/mhl.hpp>
#include <integratorxx/quadratures/treutleraldrichs.hpp>
#include <integratorxx/composite_quadratures/spherical_quadrature.hpp>
#include <gauxc/exceptions.hpp>

namespace GauXC {
namespace detail {

GridImpl::GridImpl( std::shared_ptr<quadrature_type> q, BatchSize bs ) : quad_(q) {
  generate_batcher(bs);
}

GridImpl::GridImpl( RadialQuad rq, RadialSize rs, AngularSize as, 
  RadialScale rscal, BatchSize bs )  {

  generate_stock_quadrature(rq,rs,as,rscal);
  generate_batcher(bs);

}

GridImpl::GridImpl( const GridImpl& )     = default;
GridImpl::GridImpl( GridImpl&& ) noexcept = default;

GridImpl& GridImpl::operator=( const GridImpl& )     = default;
GridImpl& GridImpl::operator=( GridImpl&& ) noexcept = default;
      
GridImpl::~GridImpl() noexcept = default;





const batcher_type& GridImpl::batcher() const { return *batcher_; }
      batcher_type& GridImpl::batcher()       { return *batcher_; }

//RadialSize  GridImpl::n_rad()        const noexcept { return n_rad_; }
//AngularSize GridImpl::n_ang()        const noexcept { return n_ang_; }
//BatchSize   GridImpl::max_batch_sz() const noexcept { return max_batch_sz_; }
//RadialScale GridImpl::rscal_factor() const noexcept { return r_scal_; }
//RadialQuad  GridImpl::radial_quad()  const noexcept { return rquad_; }


void GridImpl::generate_stock_quadrature(RadialQuad rquad, RadialSize n_rad,
  AngularSize n_ang, RadialScale r_scal) { 


  using mk_type  = IntegratorXX::MuraKnowles<double,double>;
  using mhl_type = IntegratorXX::MurrayHandyLaming<double,double>;
  using ta_type  = IntegratorXX::TreutlerAldrichs<double,double>;
  using ll_type  = IntegratorXX::LebedevLaikov<double>;

  using mk_sphere_type  = IntegratorXX::SphericalQuadrature<mk_type, ll_type>;
  using mhl_sphere_type = IntegratorXX::SphericalQuadrature<mhl_type, ll_type>;
  using ta_sphere_type  = IntegratorXX::SphericalQuadrature<ta_type, ll_type>;

  // Create Angular Quadrature
  ll_type ang_quad( n_ang.get() );

  switch( rquad ) {

    case RadialQuad::MuraKnowles:

      quad_ = std::make_shared<mk_sphere_type>(
        mk_type( n_rad.get(), r_scal.get() ),
        std::move( ang_quad )
      );
      break;

    case RadialQuad::MurrayHandyLaming:

      quad_ = std::make_shared<mhl_sphere_type>(
        mhl_type( n_rad.get(), r_scal.get() ),
        std::move( ang_quad )
      );
      break;

    case RadialQuad::TreutlerAldrichs:

      quad_ = std::make_shared<ta_sphere_type>(
        ta_type( n_rad.get(), r_scal.get() ),
        std::move( ang_quad )
      );
      break;

    default:
      GAUXC_GENERIC_EXCEPTION("Unsupported Radial Quadrature");

  }
  

}

void GridImpl::generate_batcher(BatchSize max_batch_sz) {

  batcher_ = std::make_shared< batcher_type >( 
    max_batch_sz.get(), quad_
  );

}

}
}
