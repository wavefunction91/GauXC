#pragma once

#include <gauxc/molgrid.hpp>

namespace GauXC {
namespace detail {

class MolGridImpl {

  atomic_scal_factor_map scal_factors_;
  atomic_grid_size_map   grid_sizes_;

  atomic_grid_map molgrid_;

  void generate( RadialQuad, const Molecule& );

public:

  MolGridImpl( const atomic_grid_map& );

  MolGridImpl( RadialQuad, const atomic_grid_size_map&, 
    const atomic_scal_factor_map&, const Molecule& );

  MolGridImpl( const MolGridImpl& );
  MolGridImpl( MolGridImpl&& ) noexcept;

  ~MolGridImpl() noexcept;

  size_t natoms_uniq() const;
  const Grid& get_grid( AtomicNumber ) const;
  RadialScale get_rscal_factor( AtomicNumber ) const;
  GridSize    get_grid_size( AtomicNumber ) const;
  RadialQuad  get_radial_quad( AtomicNumber ) const;

};

}
}
