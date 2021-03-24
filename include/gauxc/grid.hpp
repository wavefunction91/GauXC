#pragma once

#include <memory>
#include <gauxc/types.hpp>

namespace GauXC {

/// A named type pertaining to the size of a radial quadrature
using RadialSize   = detail::NamedType< int64_t, struct RadialSizeType  >;

/// A named type pertaining to the size of an angular quadrature
using AngularSize  = detail::NamedType< int64_t, struct AngularSizeType >;

/// A named type pertaining to the number of grid points in a quadrature batch
using BatchSize    = detail::NamedType< int64_t, struct BatchSizeType   >;

/// A named type pertaining to a scaling factor for a radial quadrature
using RadialScale  = detail::NamedType< double,  struct RadialScaleType >;

/// A type alias for the specficiation of the radial and angular dimensions of
/// a spherical quadrature
using GridSize = std::pair< RadialSize, AngularSize >;

namespace detail {
  /**
   *  @brief A class which contains the implementation details of a Grid instance
   */
  class GridImpl;
}

/**
 *  @brief A class to manage a particular spherical (atomic) quadrature
 */
class Grid {

  std::shared_ptr<detail::GridImpl> pimpl_; 
    ///< Implementation details of this particular Grid instance

public:

  Grid() = delete;

  /**
   *  @brief Construct a Grid object
   *  
   *  @param[in] rq   Radial quadrature
   *  @param[in] rs   Number of radial quadrature points
   *  @param[in] as   Number of angular quadrature points
   *  @param[in] scal Radial scaling factor ( r -> scal*r )
   *  @param[in] bs   Max number of quadrature points per batch
   */
  Grid( RadialQuad rq, RadialSize rs, AngularSize as, RadialScale scal, 
        BatchSize bs); 

  /**
   *  @brief Construct a Grid object
   *
   *  Default batch size
   *  
   *  @param[in] rq   Radial quadrature
   *  @param[in] rs   Number of radial quadrature points
   *  @param[in] as   Number of angular quadrature points
   *  @param[in] scal Radial scaling factor ( r -> scal*r )
   */
  Grid( RadialQuad rq, RadialSize rs, AngularSize as, RadialScale scal ); 

  /**
   *  @brief Construct a Grid object
   *
   *  Default radial quadrature (MuraKnowles)
   *  
   *  @param[in] rs   Number of radial quadrature points
   *  @param[in] as   Number of angular quadrature points
   *  @param[in] scal Radial scaling factor ( r -> scal*r )
   *  @param[in] bs   Max number of quadrature points per batch
   */
  Grid( RadialSize rs, AngularSize as, RadialScale scal, BatchSize bs ); 

  /**
   *  @brief Construct a Grid object
   *
   *  Default radial quadrature (MuraKnowles) and batch size
   *  
   *  @param[in] rs   Number of radial quadrature points
   *  @param[in] as   Number of angular quadrature points
   *  @param[in] scal Radial scaling factor ( r -> scal*r )
   */
  Grid( RadialSize rs, AngularSize rs, RadialScale scal ); 




  /**
   *  @brief Construct a Grid object
   *  
   *  @param[in] rq   Radial quadrature
   *  @param[in] gs   Grid dimensions
   *  @param[in] scal Radial scaling factor ( r -> scal*r )
   *  @param[in] bs   Max number of quadrature points per batch
   */
  Grid( RadialQuad rq, GridSize gs, RadialScale scal, BatchSize bs );

  /**
   *  @brief Construct a Grid object
   *
   *  Default batch size
   *  
   *  @param[in] rq   Radial quadrature
   *  @param[in] gs   Grid dimensions
   *  @param[in] scal Radial scaling factor ( r -> scal*r )
   */
  Grid( RadialQuad rq, GridSize gs, RadialScale scal );

  /**
   *  @brief Construct a Grid object
   *
   *  Default radial quadrature (MuraKnowles)
   *  
   *  @param[in] gs   Grid dimensions
   *  @param[in] scal Radial scaling factor ( r -> scal*r )
   *  @param[in] bs   Max number of quadrature points per batch
   */
  Grid( GridSize gs, RadialScale scal, BatchSize bs );

  /**
   *  @brief Construct a Grid object
   *
   *  Default radial quadrature (MuraKnowles) and batch size
   *  
   *  @param[in] gs   Grid dimensions
   *  @param[in] scal Radial scaling factor ( r -> scal*r )
   */
  Grid( GridSize gs, RadialScale scal );


  /// Copy a Grid object
  Grid( const Grid& );

  /// Move a Grid object
  Grid( Grid&& ) noexcept;

  /// Copy-assign a Grid object
  Grid& operator=( const Grid& );

  /// Move-assign a Grid object
  Grid& operator=( Grid&& ) noexcept;

  /// Destroy a Grid object
  ~Grid() noexcept;

  /**
   *  @brief Get batcher instance for underlying Grid implementation
   *
   *  Const variant
   *
   *  @returns Batcher instance pertaining to the Grid object
   */
  const batcher_type& batcher() const;

  /**
   *  @brief Get batcher instance for underlying Grid implementation
   *
   *  Non-const variant
   *
   *  @returns Batcher instance pertaining to the Grid object
   */
   batcher_type& batcher();

  
  /// Return the number of radial points in the spherical quadrature
  RadialSize  n_rad()        const noexcept;

  /// Return the number of angular points in the spherical quadrature
  AngularSize n_ang()        const noexcept;

  /// Return the max batch size for the spherical quadrature
  BatchSize   max_batch_sz() const noexcept;

  /// Return radial scaling factor for the spherical quadrature
  RadialScale rscal_factor() const noexcept; 

  /// Return the radial quadrature specifier
  RadialQuad  radial_quad()  const noexcept;

}; // class Grid

} // namespace GauXC
