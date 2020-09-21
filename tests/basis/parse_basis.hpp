#include <gauxc/basisset.hpp>
#include <gauxc/molecule.hpp>
#include <string>

namespace GauXC {

BasisSet<double> parse_basis( const Molecule& mol,
                              std::string     fname,
                              SphericalType   sph    );

}
