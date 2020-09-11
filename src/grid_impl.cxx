#include "grid_impl.hpp"

#include <integratorxx/quadratures/lebedev_laikov.hpp>
#include <integratorxx/quadratures/muraknowles.hpp>
#include <integratorxx/quadratures/mhl.hpp>
#include <integratorxx/composite_quadratures/spherical_quadrature.hpp>

namespace GauXC {
namespace detail {

GridImpl::GridImpl( RadialQuad rq, RadialSize rs, AngularSize as, 
  RadialScale rscal, BatchSize bs ) : 
  n_rad_(rs), n_ang_(as), max_batch_sz_(bs), r_scal_(rscal), rquad_(rq)  {

  generate();

}

GridImpl::GridImpl( const GridImpl& )     = default;
GridImpl::GridImpl( GridImpl&& ) noexcept = default;

GridImpl& GridImpl::operator=( const GridImpl& )     = default;
GridImpl& GridImpl::operator=( GridImpl&& ) noexcept = default;
      
GridImpl::~GridImpl() noexcept = default;





const batcher_type& GridImpl::batcher() const { return *batcher_; }
      batcher_type& GridImpl::batcher()       { return *batcher_; }

RadialSize  GridImpl::n_rad()        const noexcept { return n_rad_; }
AngularSize GridImpl::n_ang()        const noexcept { return n_ang_; }
BatchSize   GridImpl::max_batch_sz() const noexcept { return max_batch_sz_; }
RadialScale GridImpl::rscal_factor() const noexcept { return r_scal_; }
RadialQuad  GridImpl::radial_quad()  const noexcept { return rquad_; }


void GridImpl::generate() { 


  using mk_type  = IntegratorXX::MuraKnowles<double,double>;
  using mhl_type = IntegratorXX::MurrayHandyLaming<double,double>;
  using ll_type  = IntegratorXX::LebedevLaikov<double>;

  using mk_sphere_type  = IntegratorXX::SphericalQuadrature<mk_type, ll_type>;
  using mhl_sphere_type = IntegratorXX::SphericalQuadrature<mhl_type, ll_type>;

  // Create Angular Quadrature
  ll_type ang_quad( n_ang_.get() );

  switch( rquad_ ) {

    case RadialQuad::MuraKnowles:

      quad_ = std::make_shared<mk_sphere_type>(
        mk_type( n_rad_.get(), r_scal_.get() ),
        std::move( ang_quad )
      );

      batcher_ = std::make_shared< batcher_type >( 
        max_batch_sz_.get(), 
        dynamic_cast<const mk_sphere_type&>(*quad_)
      );
      break;

    case RadialQuad::MurrayHandyLaming:

      quad_ = std::make_shared<mhl_sphere_type>(
        mhl_type( n_rad_.get(), r_scal_.get() ),
        std::move( ang_quad )
      );

      batcher_ = std::make_shared< batcher_type >( 
        max_batch_sz_.get(), 
        dynamic_cast<const mhl_sphere_type&>(*quad_)
      );
      break;

    default:
      throw std::runtime_error("Unsupported Radial Quadrature");

  }
  

}

}
}
