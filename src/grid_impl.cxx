#include "grid_impl.hpp"

namespace GauXC {
namespace detail {

GridImpl::GridImpl( std::shared_ptr<quadrature_type> q, BatchSize bs ) : quad_(q) {
  generate_batcher(bs);
}

GridImpl::GridImpl( const GridImpl& )     = default;
GridImpl::GridImpl( GridImpl&& ) noexcept = default;

GridImpl& GridImpl::operator=( const GridImpl& )     = default;
GridImpl& GridImpl::operator=( GridImpl&& ) noexcept = default;
      
GridImpl::~GridImpl() noexcept = default;

const batcher_type& GridImpl::batcher() const { return *batcher_; }
      batcher_type& GridImpl::batcher()       { return *batcher_; }

void GridImpl::generate_batcher(BatchSize max_batch_sz) {

  batcher_ = std::make_shared< batcher_type >( 
    max_batch_sz.get(), quad_
  );

}

}
}
