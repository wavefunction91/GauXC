#include "catch2/catch.hpp"
#include <gauxc/molecule.hpp>
#include <gauxc/molmeta.hpp>

#include "standards.hpp"

#include <random>

using namespace GauXC;

TEST_CASE("Atom", "[moltypes]") {

  double x = 0.2, y = 0.4, z = 6.5;
  int64_t Z = 10;

  Atom atom{ AtomicNumber(Z), x, y, z };

  CHECK( atom.Z.get() == Z );
  CHECK( atom.x       == x );
  CHECK( atom.y       == y );
  CHECK( atom.z       == z );

}

TEST_CASE("Molecule", "[moltypes]") {

  std::default_random_engine gen;
  std::uniform_real_distribution<> real_dist( -10., 10 );
  std::uniform_int_distribution<>  int_dist(1,10);

  size_t natoms_gen = 40;

  SECTION("From std::vector<Atom>") {

    std::vector<Atom> atoms;
    
    for( auto i = 0; i < natoms_gen; ++i )
      atoms.push_back( {
        AtomicNumber(int_dist(gen)),
        real_dist(gen),
        real_dist(gen),
        real_dist(gen)
      });

    SECTION("Copy") {

      Molecule mol( atoms );

      CHECK( mol.natoms() == natoms_gen );
      CHECK( atoms.size() == natoms_gen );

      for( auto i = 0; i < natoms_gen; ++i ) 
        CHECK( atoms[i] == mol[i] );
      
    }

    SECTION("Move") {

      std::vector<Atom> atoms_copy(atoms);

      Molecule mol( std::move(atoms) );

      CHECK( mol.natoms() == natoms_gen );
      CHECK( atoms.size() == 0 );

      for( auto i = 0; i < natoms_gen; ++i ) 
        CHECK( atoms_copy[i] == mol[i] );
      
    }

  }

  

  SECTION("Inplace Construction") {

    Molecule mol;
    std::vector<Atom> atoms;

    for(auto i = 0ul; i < natoms_gen; ++i ) {

      AtomicNumber Z( int_dist(gen) );
      double x(real_dist(gen));
      double y(real_dist(gen));
      double z(real_dist(gen));

      atoms.push_back({Z,x,y,z});
      mol.push_back({Z,x,y,z});

    }

    CHECK( mol.natoms() == natoms_gen );
    
    for( auto i = 0; i < natoms_gen; ++i ) 
      CHECK( atoms[i] == mol[i] );

  }


}


TEST_CASE( "MolMeta", "[moltypes]" ) {

  // Water
  Molecule mol = make_water();
  MolMeta meta( mol );

  std::vector<double> rab          = meta.rab();
  std::vector<double> dist_nearest = meta.dist_nearest();

  std::vector<double> rab_ref = {
    0.00000000000, 2.68755847909, 4.34922211156,
    2.68755847909, 0.00000000000, 2.68755847909,
    4.34922211156, 2.68755847909, 0.00000000000 
  };

  for( auto i = 0; i < mol.natoms() * mol.natoms(); ++i )
    CHECK( rab[i] == Approx(rab_ref[i]) );

  for( auto i = 0; i < mol.natoms(); ++i )
    CHECK( dist_nearest[i] == Approx(2.68755847909) );
  

}
