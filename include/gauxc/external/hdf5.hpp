#include <gauxc/shell.hpp>
#include <gauxc/basisset.hpp>
#include <gauxc/molecule.hpp>


namespace GauXC {

void write_hdf5_record( const Shell<double>& shell, std::string fname, 
  std::string dset );

void write_hdf5_record( const std::vector<Shell<double>>& shell, std::string fname, 
  std::string dset );

void write_hdf5_record( const Atom& mol, std::string fname, std::string dset );
void write_hdf5_record( const Molecule& mol, std::string fname, std::string dset );
}

