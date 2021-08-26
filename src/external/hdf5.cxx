#include <gauxc/external/hdf5.hpp>
#include <highfive/H5File.hpp>
#include <stdexcept>

namespace GauXC {

using namespace HighFive;

struct shell_t {
  int32_t nprim, l, pure;
  Shell<double>::prim_array alpha, coeff;
  Shell<double>::cart_array O;
};

hid_t create_shell_type() {

  hsize_t prim_dims[1] = {16};
  hsize_t cart_dims[1] = {3};
  auto prim_array_type = H5Tarray_create( H5T_NATIVE_DOUBLE, 1, prim_dims );
  auto cart_array_type = H5Tarray_create( H5T_NATIVE_DOUBLE, 1, cart_dims );

  auto shell_type = H5Tcreate( H5T_COMPOUND, sizeof(shell_t) );
  H5Tinsert( shell_type, "NPRIM",  HOFFSET( shell_t, nprim ), H5T_NATIVE_INT );
  H5Tinsert( shell_type, "L",      HOFFSET( shell_t, l ),     H5T_NATIVE_INT );
  H5Tinsert( shell_type, "PURE",   HOFFSET( shell_t, pure ),  H5T_NATIVE_INT );
  H5Tinsert( shell_type, "ALPHA",  HOFFSET( shell_t, alpha),  prim_array_type );
  H5Tinsert( shell_type, "COEFF",  HOFFSET( shell_t, coeff),  prim_array_type );
  H5Tinsert( shell_type, "ORIGIN", HOFFSET( shell_t, O),      cart_array_type );
  
  return shell_type;
}

void write_hdf5_record( const Shell<double>& shell, std::string fname, 
  std::string dset ) {

  File file( fname, File::OpenOrCreate );
  
  auto shell_type = create_shell_type();

  DataSpace space(1);
  auto d_id = H5Dcreate( file.getId(), dset.c_str(), shell_type, space.getId(),
    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

  if( d_id < 0 ) throw std::runtime_error("Dataset Creation Failed");

  shell_t sh_ {
    shell.nprim(), shell.l(), shell.pure(),
    shell.alpha(), shell.coeff(), shell.O()
  };

  H5Dwrite( d_id, shell_type, space.getId(), space.getId(), H5P_DEFAULT, 
    &sh_ );

  H5Tclose( shell_type );
  H5Dclose( d_id );

}



void write_hdf5_record( const std::vector<Shell<double>>& basis, std::string fname, 
  std::string dset ) {


  File file( fname, File::OpenOrCreate );
  
  auto shell_type = create_shell_type();

  DataSpace space(basis.size());
  auto d_id = H5Dcreate( file.getId(), dset.c_str(), shell_type, space.getId(),
    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

  if( d_id < 0 ) throw std::runtime_error("Dataset Creation Failed");

  std::vector<shell_t> shells;
  for( auto& shell : basis ) {
    shells.push_back(
      shell_t{
        shell.nprim(), shell.l(), shell.pure(),
        shell.alpha(), shell.coeff(), shell.O()
      });
  }

  H5Dwrite( d_id, shell_type, space.getId(), space.getId(), H5P_DEFAULT, 
    shells.data() );

  H5Tclose( shell_type );
  H5Dclose( d_id );


}

hid_t create_atom_type() {

  auto atom_type = H5Tcreate( H5T_COMPOUND, sizeof(Atom) );
  H5Tinsert( atom_type, "Atomic Number", HOFFSET( Atom, Z), H5T_NATIVE_INT    );
  H5Tinsert( atom_type, "X Coordinate", HOFFSET( Atom, x ), H5T_NATIVE_DOUBLE );
  H5Tinsert( atom_type, "Y Coordinate", HOFFSET( Atom, y ), H5T_NATIVE_DOUBLE );
  H5Tinsert( atom_type, "Z Coordinate", HOFFSET( Atom, z ), H5T_NATIVE_DOUBLE );

  return atom_type;

}



void write_hdf5_record( const Atom& atom, std::string fname, std::string dset ) {

  File file( fname, File::OpenOrCreate );
  
  auto atom_type = create_atom_type();

  DataSpace space(1);
  auto d_id = H5Dcreate( file.getId(), dset.c_str(), atom_type, space.getId(),
    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

  if( d_id < 0 ) throw std::runtime_error("Dataset Creation Failed");

  H5Dwrite( d_id, atom_type, space.getId(), space.getId(), H5P_DEFAULT, 
    &atom );

  H5Tclose( atom_type );
  H5Dclose( d_id );

}

void write_hdf5_record( const Molecule& mol, std::string fname, std::string dset ) {

  File file( fname, File::OpenOrCreate );
  
  auto atom_type = create_atom_type();

  DataSpace space(mol.size());
  auto d_id = H5Dcreate( file.getId(), dset.c_str(), atom_type, space.getId(),
    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

  if( d_id < 0 ) throw std::runtime_error("Dataset Creation Failed");

  H5Dwrite( d_id, atom_type, space.getId(), space.getId(), H5P_DEFAULT, 
    mol.data() );

  H5Tclose( atom_type );
  H5Dclose( d_id );

}


}
