#pragma once

#include <memory>
#include <gauxc/types.hpp>

namespace GauXC {


enum class RadialQuad {
  MuraKnowles,
  MurrayHandyLaming
};


using RadialSize   = detail::NamedType< int64_t, struct RadialSizeType  >;
using AngularSize  = detail::NamedType< int64_t, struct AngularSizeType >;
using BatchSize    = detail::NamedType< int64_t, struct BatchSizeType   >;
using RadialScale  = detail::NamedType< double,  struct RadialScaleType >;

using GridSize = std::pair< RadialSize, AngularSize >;

namespace detail {
  class GridImpl;
}


class Grid {

  std::shared_ptr<detail::GridImpl> pimpl_;

public:

  Grid() = delete;

  Grid( RadialQuad, RadialSize, AngularSize, RadialScale, BatchSize); 
  Grid( RadialQuad, RadialSize, AngularSize, RadialScale ); 
  Grid( RadialSize, AngularSize, RadialScale, BatchSize ); 
  Grid( RadialSize, AngularSize, RadialScale ); 

  Grid( RadialQuad, GridSize, RadialScale, BatchSize );
  Grid( RadialQuad, GridSize, RadialScale );
  Grid( GridSize, RadialScale, BatchSize );
  Grid( GridSize, RadialScale );


  Grid( const Grid& );
  Grid( Grid&& ) noexcept;

  Grid& operator=( const Grid& );
  Grid& operator=( Grid&& ) noexcept;

  ~Grid() noexcept;

  const batcher_type& batcher() const;
        batcher_type& batcher()      ;

  RadialSize  n_rad()        const noexcept;
  AngularSize n_ang()        const noexcept;
  BatchSize   max_batch_sz() const noexcept;
  RadialScale rscal_factor() const noexcept; 

  RadialQuad  radial_quad()  const noexcept;

};

}
