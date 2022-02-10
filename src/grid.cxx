#include "grid_impl.hpp"


namespace GauXC {

constexpr BatchSize default_batch_sz = BatchSize(512);

Grid::Grid( std::shared_ptr<quadrature_type> q ) :
  pimpl_( std::make_shared<detail::GridImpl>(q,default_batch_sz) ) { }

Grid::Grid( const Grid& )     = default;
Grid::Grid( Grid&& ) noexcept = default;

Grid& Grid::operator=( const Grid& )     = default;
Grid& Grid::operator=( Grid&& ) noexcept = default;
      
Grid::~Grid() noexcept = default;

const batcher_type& Grid::batcher() const { return pimpl_->batcher(); }
      batcher_type& Grid::batcher()       { return pimpl_->batcher(); }

}
