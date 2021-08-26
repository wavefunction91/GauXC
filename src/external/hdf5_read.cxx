#include "hdf5_util.hpp"

namespace GauXC {

using namespace HighFive;


void read_hdf5_record( std::vector<Shell<double>>& basis, std::string fname, 
  std::string dset ) {


  File file( fname, File::ReadOnly );
  

  auto d_id = H5Dopen( file.getId(), dset.c_str(), H5P_DEFAULT );
  if( d_id < 0 ) throw std::runtime_error("Dataset Open Failed");

  auto space_id = H5Dget_space( d_id );
  if( space_id < 0 ) throw std::runtime_error( "Space Retreival failed" );

  auto ndims = H5Sget_simple_extent_ndims( space_id );
  if( ndims != 1 ) throw std::runtime_error("Only supported for 1D data structures");

  hsize_t size;
  H5Sget_simple_extent_dims( space_id, &size, NULL );


  std::vector<shell_t> shells( size );
  auto shell_type = create_shell_type();
  H5Dread( d_id, shell_type, space_id, space_id, H5P_DEFAULT, shells.data() );

  basis.resize( size );
  for( auto i = 0ul; i < size; ++i ) {
    auto& sh = shells[i];
    basis[i] = Shell<double>( PrimSize(sh.nprim), AngularMomentum(sh.l),
      SphericalType(sh.pure), sh.alpha, sh.coeff, sh.O, false );
  }



  H5Tclose( shell_type );
  H5Dclose( d_id );
  H5Sclose( space_id );


}




void read_hdf5_record( std::vector<Atom>& mol, std::string fname, std::string dset ) {

  File file( fname, File::ReadOnly );
  
  auto d_id = H5Dopen( file.getId(), dset.c_str(), H5P_DEFAULT );
  if( d_id < 0 ) throw std::runtime_error("Dataset Open Failed");

  auto space_id = H5Dget_space( d_id );
  if( space_id < 0 ) throw std::runtime_error( "Space Retreival failed" );

  auto ndims = H5Sget_simple_extent_ndims( space_id );
  if( ndims != 1 ) throw std::runtime_error("Only supported for 1D data structures");

  hsize_t size;
  H5Sget_simple_extent_dims( space_id, &size, NULL );


  auto atom_type = create_atom_type();
  mol.resize(size);
  H5Dread( d_id, atom_type, space_id, space_id, H5P_DEFAULT, mol.data() );

  H5Tclose( atom_type );
  H5Dclose( d_id );
  H5Sclose( space_id );
}


}
