#include <gauxc/gauxc_config.hpp>
#ifdef GAUXC_ENABLE_HDF5
#include <gauxc/shell.hpp>
#include <gauxc/atom.hpp>

namespace GauXC {
void write_hdf5_record( const std::vector<Shell<double>>& shell, std::string fname, std::string dset );
void write_hdf5_record( const std::vector<Atom>& mol, std::string fname, std::string dset );
void read_hdf5_record( std::vector<Shell<double>>& shell, std::string fname, std::string dset );
void read_hdf5_record( std::vector<Atom>& mol, std::string fname, std::string dset );
}
#endif
