#include "catch2/catch.hpp"
#include <gauxc/basisset.hpp>
#include <gauxc/molecule.hpp>

#include <random>

using namespace GauXC;

//H     0
//S   3   1.00
//      0.1873113696D+02       0.3349460434D-01
//      0.2825394365D+01       0.2347269535D+00
//      0.6401216923D+00       0.8137573261D+00
//S   1   1.00
//      0.1612777588D+00       1.0000000
//****
//O     0
//S   6   1.00
//      0.5484671660D+04       0.1831074430D-02
//      0.8252349460D+03       0.1395017220D-01
//      0.1880469580D+03       0.6844507810D-01
//      0.5296450000D+02       0.2327143360D+00
//      0.1689757040D+02       0.4701928980D+00
//      0.5799635340D+01       0.3585208530D+00
//SP   3   1.00
//      0.1553961625D+02      -0.1107775495D+00       0.7087426823D-01
//      0.3599933586D+01      -0.1480262627D+00       0.3397528391D+00
//      0.1013761750D+01       0.1130767015D+01       0.7271585773D+00
//SP   1   1.00
//      0.2700058226D+00       0.1000000000D+01       0.1000000000D+01
//D   1   1.00
//      0.8000000000D+00       1.0000000

TEST_CASE("BasisSet", "[basisset]") {

  BasisSet<double> basis;

  Molecule mol;
  mol.emplace_back(AtomicNumber(1), 0., 1.579252144093028,  2.174611055780858);
  mol.emplace_back(AtomicNumber(8), 0., 0.000000000000000,  0.000000000000000);
  mol.emplace_back(AtomicNumber(1), 0., 1.579252144093028, -2.174611055780858);

  bool test_spherical = false;
  SECTION( "Cartesian" ) { test_spherical = false; }
  SECTION( "Spherical" ) { test_spherical = true;  }

  for( const auto& atom : mol ) {

    if( atom.Z == 1ll ) {

      basis.emplace_back(
        Shell<double>( 
          PrimSize(3), AngularMomentum(0), SphericalType(false),
          {0.3349460434e-01, 0.2347269535e+00, 0.8137573261e+00},
          {0.1873113696e+02, 0.2825394365e+01, 0.6401216923e+00},
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>( 
          PrimSize(1), AngularMomentum(0), SphericalType(false),
          {0.1612777588},
          {1.0000000000},
          { atom.x, atom.y, atom.z }
        )
      );

    } else {

      basis.emplace_back(
        Shell<double>( 
          PrimSize(6), AngularMomentum(0), SphericalType(false),
          { 0.5484671660e+04, 0.8252349460e+03, 0.1880469580e+03, 
            0.5296450000e+02, 0.1689757040e+02, 0.5799635340e+01 }, 
          { 0.1831074430e-02, 0.1395017220e-01, 0.6844507810e-01,
            0.2327143360e+00, 0.4701928980e+00, 0.3585208530e+00 },
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>( 
          PrimSize(3), AngularMomentum(0), SphericalType(false),
          { 0.1553961625e+02,  0.3599933586e+01, 0.1013761750e+01},
          {-0.1107775495e+00, -0.1480262627e+00, 0.1130767015e+01},
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>( 
          PrimSize(3), AngularMomentum(1), SphericalType(false),
          { 0.1553961625e+02,  0.3599933586e+01, 0.1013761750e+01},
          { 0.7087426823e-01,  0.3397528391e+00, 0.7271585773e+00},
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>( 
          PrimSize(1), AngularMomentum(0), SphericalType(false),
          { 0.2700058226 },
          { 1.0000000000 },
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>( 
          PrimSize(1), AngularMomentum(1), SphericalType(false),
          { 0.2700058226 },
          { 1.0000000000 },
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>( 
          PrimSize(1), AngularMomentum(2), SphericalType(test_spherical),
          { 0.8 },
          { 1.0 },
          { atom.x, atom.y, atom.z }
        )
      );

    }

  }

  CHECK( basis.nshells() == 10 );
  CHECK( basis.nbf()     == (test_spherical ? 18 : 19) );
  basis.generate_shell_to_ao();


  std::vector<int32_t> ref_shell_to_ao = {
  0, 1, // H1
  2, 3, 4, 7, 8, 11, // O
  (test_spherical ? 16 : 17), (test_spherical ? 17 : 18)  // H2
  };

  CHECK( basis.shell_to_first_ao() == ref_shell_to_ao );

  for(auto i = 0; i < basis.nshells(); ++i) {
    auto [sh_st,sh_en] = basis.shell_to_ao_range(i);
    CHECK(sh_st == ref_shell_to_ao[i]);
    CHECK(sh_en == ref_shell_to_ao[i] + basis[i].size());
  }

}
