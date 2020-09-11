#include "grid_impl.hpp"


namespace GauXC {

constexpr BatchSize default_batch_sz = BatchSize(512);

Grid::Grid( RadialQuad rq, RadialSize rs, AngularSize as, RadialScale rscal,
  BatchSize bs ) :
  pimpl_(std::make_shared<detail::GridImpl>(rq, rs, as, rscal, bs)) { }

Grid::Grid( RadialQuad rq, RadialSize rs, AngularSize as, RadialScale rscal ) :
  Grid( rq, rs, as, rscal, default_batch_sz ) { }

Grid::Grid( RadialSize rs, AngularSize as, RadialScale rscal, BatchSize bs ) :
  Grid( RadialQuad::MuraKnowles, rs, as, rscal, bs ) { }

Grid::Grid( RadialSize rs, AngularSize as, RadialScale rscal ) :
  Grid( rs, as, rscal, default_batch_sz ) { }


Grid::Grid( RadialQuad rq, GridSize gs, RadialScale rscal, BatchSize bs ) :
  Grid( rq, gs.first, gs.second, rscal, bs ) { }
Grid::Grid( RadialQuad rq, GridSize gs, RadialScale rscal ) :
  Grid( rq, gs.first, gs.second, rscal ) { }
Grid::Grid( GridSize gs, RadialScale rscal, BatchSize bs ) :
  Grid( gs.first, gs.second, rscal, bs ) { }
Grid::Grid( GridSize gs, RadialScale rscal ) :
  Grid( gs.first, gs.second, rscal ) { }


Grid::Grid( const Grid& )     = default;
Grid::Grid( Grid&& ) noexcept = default;

Grid& Grid::operator=( const Grid& )     = default;
Grid& Grid::operator=( Grid&& ) noexcept = default;
      
Grid::~Grid() noexcept = default;


const batcher_type& Grid::batcher() const { return pimpl_->batcher(); }
      batcher_type& Grid::batcher()       { return pimpl_->batcher(); }

RadialSize  Grid::n_rad()        const noexcept { return pimpl_->n_rad()       ; }
AngularSize Grid::n_ang()        const noexcept { return pimpl_->n_ang()       ; }
BatchSize   Grid::max_batch_sz() const noexcept { return pimpl_->max_batch_sz(); }
RadialScale Grid::rscal_factor() const noexcept { return pimpl_->rscal_factor(); }
RadialQuad  Grid::radial_quad()  const noexcept { return pimpl_->radial_quad() ; }

}
