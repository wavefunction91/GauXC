#pragma once

#include <gauxc/grid.hpp>

namespace GauXC {
namespace detail {

class GridImpl {

  RadialSize   n_rad_;
  AngularSize  n_ang_;
  BatchSize    max_batch_sz_;
  RadialScale  r_scal_;

  RadialQuad  rquad_;

  std::shared_ptr< quadrature_type > quad_    = nullptr;
  std::shared_ptr< batcher_type    > batcher_ = nullptr;

  void generate();

public:

  GridImpl() = delete;

  GridImpl( RadialQuad, RadialSize, AngularSize, RadialScale, BatchSize ); 

  GridImpl( const GridImpl& );
  GridImpl( GridImpl&& ) noexcept;

  GridImpl& operator=( const GridImpl& );
  GridImpl& operator=( GridImpl&& ) noexcept;

  ~GridImpl() noexcept;
  
  const batcher_type& batcher() const;
        batcher_type& batcher()      ;

  RadialSize  n_rad()        const noexcept;
  AngularSize n_ang()        const noexcept;
  BatchSize   max_batch_sz() const noexcept;
  RadialScale rscal_factor() const noexcept; 

  RadialQuad  radial_quad()  const noexcept;

};

}
}
