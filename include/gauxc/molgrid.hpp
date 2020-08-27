#pragma once

#include <gauxc/molecule.hpp>
#include <gauxc/grid.hpp>

#include <unordered_map>

namespace GauXC {

using atomic_scal_factor_map =
  std::unordered_map< AtomicNumber, RadialScale >;

using atomic_grid_size_map =
  std::unordered_map< AtomicNumber, GridSize >;

using atomic_grid_map = std::unordered_map< AtomicNumber, Grid >;

namespace detail {
  class MolGridImpl;
}

class MolGrid {

  std::shared_ptr<detail::MolGridImpl> pimpl_;

public:

  MolGrid( const atomic_grid_map& );

  MolGrid( RadialQuad, const atomic_grid_size_map&, const atomic_scal_factor_map&, 
    const Molecule& );
  MolGrid( RadialQuad, const atomic_grid_size_map&, const Molecule& );

  MolGrid( RadialQuad, AtomicGridSizeDefault, const atomic_scal_factor_map&, 
    const Molecule& );
  MolGrid( RadialQuad, AtomicGridSizeDefault, const Molecule& );

  MolGrid( const atomic_grid_size_map&, const atomic_scal_factor_map&, 
    const Molecule& );
  MolGrid( const atomic_grid_size_map&, const Molecule& );

  MolGrid( AtomicGridSizeDefault, const atomic_scal_factor_map&, 
    const Molecule& );
  MolGrid( AtomicGridSizeDefault, const Molecule& );

  MolGrid( const MolGrid& );
  MolGrid( MolGrid&& ) noexcept;

  ~MolGrid() noexcept;

  size_t natoms_uniq() const;
  const Grid& get_grid( AtomicNumber ) const;
        Grid& get_grid( AtomicNumber )      ;
  RadialScale get_rscal_factor( AtomicNumber ) const;
  GridSize    get_grid_size( AtomicNumber ) const ;
  RadialQuad  get_radial_quad( AtomicNumber ) const;

  size_t max_nbatches() const;

};

}
