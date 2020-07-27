#include <gauxc/basisset.hpp>
#include <gauxc/molecule.hpp>

namespace GauXC {

Molecule         make_water();
Molecule         make_benzene();
BasisSet<double> make_631Gd( const Molecule&, SphericalType );
BasisSet<double> make_ccpvdz( const Molecule&, SphericalType );

}
