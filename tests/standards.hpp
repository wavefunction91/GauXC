#include <gauxc/basisset.hpp>
#include <gauxc/molecule.hpp>

namespace GauXC {

Molecule         make_water();
BasisSet<double> make_631Gd( const Molecule&, SphericalType );

}
