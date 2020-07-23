#include "molgrid_impl.hpp"
#include <gauxc/molgrid/defaults.hpp>

namespace GauXC {

MolGrid::MolGrid( const atomic_grid_map& ag ) :
  pimpl_( std::make_shared<detail::MolGridImpl>( ag ) ) { }
  
MolGrid::MolGrid( 
  RadialQuad                    rq,
  const atomic_grid_size_map&   grid_sz, 
  const atomic_scal_factor_map& rad_scl, 
  const Molecule& mol 
) : pimpl_( std::make_shared<detail::MolGridImpl>( rq, grid_sz, rad_scl, mol ) ) { }

MolGrid::MolGrid( 
  RadialQuad                    rq,
  const atomic_grid_size_map&   grid_sz, 
  const Molecule& mol 
) : MolGrid( rq, grid_sz, get_default_scaling_factors(rq, mol.maxZ()), mol ) { }



MolGrid::MolGrid( 
  RadialQuad                    rq,
  AtomicGridSizeDefault         grid_sz, 
  const atomic_scal_factor_map& rad_scl, 
  const Molecule& mol 
) : MolGrid( rq, get_default_grid_sizes(grid_sz, mol.maxZ()), rad_scl, mol ) { }

MolGrid::MolGrid( 
  RadialQuad                    rq,
  AtomicGridSizeDefault         grid_sz, 
  const Molecule& mol 
) : MolGrid( rq, get_default_grid_sizes(grid_sz, mol.maxZ()), mol ) { }


MolGrid::MolGrid( 
  const atomic_grid_size_map&   grid_sz, 
  const atomic_scal_factor_map& rad_scl, 
  const Molecule& mol
) : MolGrid( RadialQuad::MuraKnowles, grid_sz, rad_scl, mol ) { } 

MolGrid::MolGrid( 
  const atomic_grid_size_map&   grid_sz, 
  const Molecule& mol
) : MolGrid( RadialQuad::MuraKnowles, grid_sz, mol ) { } 

MolGrid::MolGrid( 
  AtomicGridSizeDefault   grid_sz, 
  const atomic_scal_factor_map& rad_scl, 
  const Molecule& mol
) : MolGrid( RadialQuad::MuraKnowles, grid_sz, rad_scl, mol ) { } 

MolGrid::MolGrid( 
  AtomicGridSizeDefault   grid_sz, 
  const Molecule& mol
) : MolGrid( RadialQuad::MuraKnowles, grid_sz, mol ) { } 

MolGrid::MolGrid( const MolGrid& )     = default;
MolGrid::MolGrid( MolGrid&& ) noexcept = default;
MolGrid::~MolGrid() noexcept = default;

size_t MolGrid::natoms_uniq() const { return pimpl_->natoms_uniq(); }
const Grid& MolGrid::get_grid( AtomicNumber Z ) const { 
  return pimpl_->get_grid(Z); 
}
RadialScale MolGrid::get_rscal_factor( AtomicNumber Z ) const {
  return pimpl_->get_rscal_factor(Z);
}
GridSize MolGrid::get_grid_size( AtomicNumber Z ) const {
  return pimpl_->get_grid_size(Z);
}
RadialQuad MolGrid::get_radial_quad( AtomicNumber Z ) const {
  return pimpl_->get_radial_quad(Z);
}

}

