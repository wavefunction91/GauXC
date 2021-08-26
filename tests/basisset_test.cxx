#include "ut_common.hpp"
#include "catch2/catch.hpp"
#include <gauxc/basisset.hpp>
#include <gauxc/basisset_map.hpp>
#include <gauxc/molecule.hpp>
#include <gauxc/external/hdf5.hpp>

#include "standards.hpp"

#include <random>
#include <algorithm>

#include <gauxc/gauxc_config.hpp>

#ifdef GAUXC_ENABLE_MPI
#include <mpi.h>
#endif

using namespace GauXC;



TEST_CASE("Shell", "[basisset]") {

  using prim_array = Shell<double>::prim_array;
  using cart_array = Shell<double>::cart_array;

  const cart_array center = {0., 1., 0.};

  auto cutFunc = []( double alpha, double tol ) -> double {
    const double log_tol  = -std::log(tol);
    const double log_alph =  std::log(alpha);
    return std::sqrt( (log_tol + 0.5 * log_alph)/alpha );
  };


  const double sqrt_pi = std::sqrt(M_PI);
  auto s_int = [=](double a) { return sqrt_pi / std::sqrt(a); };
  auto p_int = [=](double a) { return sqrt_pi / (2.*std::pow(a,1.5)); };
  auto d_int = [=](double a) { return 3.*sqrt_pi / (4.*std::pow(a,2.5)); };

  SECTION("Single Gaussian") {

    const auto nprim = PrimSize(1);
    const prim_array alpha = {0.8};
    const prim_array coeff = {0.5};

    SECTION("S Function") {

      Shell<double> sh( nprim, AngularMomentum(0), SphericalType(false),
        alpha, coeff, center );

      const double ncoeff = 1./std::sqrt(std::pow(s_int(2*alpha[0]),3.));

      CHECK( sh.nprim() == nprim );
      CHECK( sh.l()     == 0     );
      CHECK( !sh.pure()  );
      CHECK( sh.alpha()[0] == alpha[0] );
      CHECK( sh.coeff()[0] == Approx(ncoeff) );
      CHECK( sh.size() == 1 );


      CHECK( sh.cutoff_radius() == Approx(cutFunc(alpha[0], 1e-10)) );

      double exact_int = 0.;
      for( int32_t i = 0; i < 1; ++i )
      for( int32_t j = 0; j < 1; ++j )
        exact_int += sh.coeff()[i] * sh.coeff()[j] * 
          std::pow( s_int(sh.alpha()[i] + sh.alpha()[j]), 3. );

      CHECK( exact_int == Approx(1.) );
    }

    SECTION("P Function") {

      Shell<double> sh( nprim, AngularMomentum(1), SphericalType(false),
        alpha, coeff, center );

      const double exact_int = std::pow(s_int(2*alpha[0]),2.) *  p_int(2*alpha[0]);
      const double ncoeff = 1./std::sqrt(exact_int);

      CHECK( sh.nprim() == nprim );
      CHECK( sh.l()     == 1     );
      CHECK( !sh.pure() );
      CHECK( sh.alpha()[0] == alpha[0] );
      CHECK( sh.coeff()[0] == Approx(ncoeff) );
      CHECK( sh.size() == 3 );


      CHECK( sh.cutoff_radius() == Approx(cutFunc(alpha[0], 1e-10)) );

    }

    SECTION("D Function") {

      Shell<double> sh( nprim, AngularMomentum(2), SphericalType(false),
        alpha, coeff, center );

      const double exact_int = std::pow(s_int(2*alpha[0]),2.) *  d_int(2*alpha[0]);
      const double ncoeff = 1./std::sqrt(exact_int);

      CHECK( sh.nprim() == nprim );
      CHECK( sh.l()     == 2     );
      CHECK( !sh.pure() );
      CHECK( sh.alpha()[0] == alpha[0] );
      CHECK( sh.coeff()[0] == Approx(ncoeff) );
      CHECK( sh.size() == 6 );


      CHECK( sh.cutoff_radius() == Approx(cutFunc(alpha[0], 1e-10)) );

    }

  }

  SECTION("Multiple Gaussians") {

    const auto nprim = PrimSize(3);
    const prim_array coeff = {0.3349460434e-01, 0.2347269535e+00, 0.8137573261e+00};
    const prim_array alpha = {0.1873113696e+02, 0.2825394365e+01, 0.6401216923e+00};

    Shell<double> sh( nprim, AngularMomentum(0), SphericalType(false),
      alpha, coeff, center );

    double exact_int = 0.;
    for( int32_t i = 0; i < 3; ++i )
    for( int32_t j = 0; j < 3; ++j )
      exact_int += sh.coeff()[i] * sh.coeff()[j] * 
        std::pow( s_int(sh.alpha()[i] + sh.alpha()[j]), 3. );


    CHECK( exact_int == Approx(1.) );


    std::vector<double> indiv_cutoffs( 3 );
    std::transform( sh.alpha().begin(), sh.alpha().begin()+3, indiv_cutoffs.begin(),
      [&](double a) { return cutFunc(a, 1e-10); });

    CHECK( sh.cutoff_radius() == Approx(*std::max_element(
      indiv_cutoffs.begin(), indiv_cutoffs.end()
    )));

  }

  SECTION( "Cutoff Nondefault Tolerance" ) {

    const auto nprim = PrimSize(1);
    const prim_array alpha = {0.8};
    const prim_array coeff = {0.5};

    Shell<double> sh( nprim, AngularMomentum(2), SphericalType(false),
      alpha, coeff, center );

    CHECK(sh.cutoff_radius() == Approx( cutFunc( alpha[0], 1e-10 ) ));    

    sh.set_shell_tolerance( 1e-10 );
    CHECK(sh.cutoff_radius() == Approx( cutFunc( alpha[0], 1e-10 ) ));    

    sh.set_shell_tolerance( 1e-7 );
    CHECK(sh.cutoff_radius() == Approx( cutFunc( alpha[0], 1e-7 ) ));    

  }


  SECTION("Spherical") {

    const auto nprim = PrimSize(1);
    const prim_array alpha = {0.8};
    const prim_array coeff = {0.5};

    Shell<double> sh( nprim, AngularMomentum(2), SphericalType(true),
      alpha, coeff, center );

    CHECK( sh.size() == 5 );

  }

}


TEST_CASE("BasisSet", "[basisset]") {


  bool test_spherical = false;
  SECTION( "Cartesian" ) { test_spherical = false; }
  SECTION( "Spherical" ) { test_spherical = true;  }


  Molecule mol = make_water();
  BasisSet<double> basis = make_631Gd(mol, SphericalType(test_spherical));

  CHECK( basis.nshells() == 10 );
  CHECK( basis.nbf()     == (test_spherical ? 18 : 19) );
  BasisSetMap basis_map( basis, mol );

  std::vector<int32_t> ref_shell_to_ao = {
  0, 1, // H1
  2, 3, 4, 7, 8, 11, // O
  (test_spherical ? 16 : 17), (test_spherical ? 17 : 18)  // H2
  };

  CHECK( basis_map.shell_to_first_ao() == ref_shell_to_ao );
  CHECK( std::none_of( basis_map.shell_to_center().begin(), basis_map.shell_to_center().end(), [](auto i){ return i == -1;} ) );

  for(auto i = 0; i < basis.nshells(); ++i) {
    auto [sh_st,sh_en] = basis_map.shell_to_ao_range(i);
    CHECK(sh_st == ref_shell_to_ao[i]);
    CHECK(sh_en == ref_shell_to_ao[i] + basis[i].size());
  }

}



TEST_CASE("HDF5-BASISSET", "[basisset]") {

#ifdef GAUXC_ENABLE_MPI
  int world_rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &world_rank );
  if( world_rank ) return; // Only run on root rank
#endif


  Molecule mol = make_water();
  BasisSet<double> basis = make_631Gd(mol, SphericalType(false));
  
  // Write file
  const std::string fname = GAUXC_REF_DATA_PATH "/test_basis.hdf5";
  write_hdf5_record( basis, fname , "/BASIS" );

  // Read File
  BasisSet<double> basis_read;
  read_hdf5_record( basis_read, fname, "/BASIS" );

  // Check that IO was correct
  CHECK( basis == basis_read );

  std::remove( fname.c_str() ); // Delete the test file

}
