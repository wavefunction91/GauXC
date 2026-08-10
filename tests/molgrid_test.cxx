/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy).
 *
 * (c) 2024-2025, Microsoft Corporation
 *
 * All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include "catch2/catch.hpp"
#include <gauxc/exceptions.hpp>
#include <gauxc/molgrid.hpp>
#include <gauxc/molgrid/defaults.hpp>


#include <random>

using namespace GauXC;

TEST_CASE("MolGrid Defaults", "[molgrid]") {

  SECTION("MK Defaults") {

#if 0
    auto default_scaling = 
      get_default_scaling_factors( RadialQuad::MuraKnowles, AtomicNumber(100) );
     
    for( auto [Z, alpha] : default_scaling ) {

      switch( Z.get() ) {
        case 3:
        case 4:
        case 11:
        case 12:
        case 19:
        case 20:
        case 37:
        case 38:
        case 55:
        case 56:
        case 87:
        case 88:
          CHECK( alpha == 7.0 );
          break;
        default:
          CHECK( alpha == 5.0 );
          break;

      }

    }
#else
    for( auto i = 0; i < 100; ++i ) {
      AtomicNumber Z(i);
      auto alpha = 
        default_radial_scaling_factor(RadialQuad::MuraKnowles,Z).get();
      switch(i) {
        case 3:
        case 4:
        case 11:
        case 12:
        case 19:
        case 20:
        case 37:
        case 38:
        case 55:
        case 56:
        case 87:
        case 88:
          CHECK( alpha == 7.0 );
          break;
        default:
          CHECK( alpha == 5.0 );
          break;
      }
    }
#endif

  }


  SECTION("MHL Defaults") {


    for( auto i = 0; i < 100; ++i ) {
      AtomicNumber Z(i);
      auto alpha = 
        default_radial_scaling_factor(RadialQuad::MurrayHandyLaming,Z).get();
      switch(i) {
        case    1: CHECK( alpha == Approx(4.7243153124839748e-01)); break;
        case    2: CHECK( alpha == Approx(2.9290754937400643e-01)); break;
        case    3: CHECK( alpha == Approx(1.3700514406203526e+00)); break;
        case    4: CHECK( alpha == Approx(9.9210621562163470e-01)); break;
        case    5: CHECK( alpha == Approx(8.0313360312227566e-01)); break;
        case    6: CHECK( alpha == Approx(6.6140414374775647e-01)); break;
        case    7: CHECK( alpha == Approx(6.1416099062291674e-01)); break;
        case    8: CHECK( alpha == Approx(5.6691783749807700e-01)); break;
        case    9: CHECK( alpha == Approx(4.7243153124839748e-01)); break;
        case   10: CHECK( alpha == Approx(3.5904796374878206e-01)); break;
        case   11: CHECK( alpha == Approx(1.7007535124942308e+00)); break;
        case   12: CHECK( alpha == Approx(1.4172945937451924e+00)); break;
        case   13: CHECK( alpha == Approx(1.1810788281209936e+00)); break;
        case   14: CHECK( alpha == Approx(1.0393493687464743e+00)); break;
        case   15: CHECK( alpha == Approx(9.4486306249679497e-01)); break;
        case   16: CHECK( alpha == Approx(9.4486306249679497e-01)); break;
        case   17: CHECK( alpha == Approx(9.4486306249679497e-01)); break;
        case   18: CHECK( alpha == Approx(6.7085277437272439e-01)); break;
        case   19: CHECK( alpha == Approx(2.0786987374929486e+00)); break;
        case   20: CHECK( alpha == Approx(1.7007535124942308e+00)); break;
        case   21: CHECK( alpha == Approx(1.5117808999948719e+00)); break;
        case   22: CHECK( alpha == Approx(1.3228082874955129e+00)); break;
        case   23: CHECK( alpha == Approx(1.2755651343706731e+00)); break;
        case   24: CHECK( alpha == Approx(1.3228082874955129e+00)); break;
        case   25: CHECK( alpha == Approx(1.3228082874955129e+00)); break;
        case   26: CHECK( alpha == Approx(1.3228082874955129e+00)); break;
        case   27: CHECK( alpha == Approx(1.2755651343706731e+00)); break;
        case   28: CHECK( alpha == Approx(1.2755651343706731e+00)); break;
        case   29: CHECK( alpha == Approx(1.2755651343706731e+00)); break;
        case   30: CHECK( alpha == Approx(1.2755651343706731e+00)); break;
        case   31: CHECK( alpha == Approx(1.2283219812458335e+00)); break;
        case   32: CHECK( alpha == Approx(1.1810788281209936e+00)); break;
        case   33: CHECK( alpha == Approx(1.0865925218713142e+00)); break;
        case   34: CHECK( alpha == Approx(1.0865925218713142e+00)); break;
        case   35: CHECK( alpha == Approx(1.0865925218713142e+00)); break;
        case   36: CHECK( alpha == Approx(8.3147949499717955e-01)); break;
        case   37: CHECK( alpha == Approx(2.2204281968674682e+00)); break;
        case   38: CHECK( alpha == Approx(1.8897261249935899e+00)); break;
        case   39: CHECK( alpha == Approx(1.7007535124942308e+00)); break;
        case   40: CHECK( alpha == Approx(1.4645377468700322e+00)); break;
        case   41: CHECK( alpha == Approx(1.3700514406203526e+00)); break;
        case   42: CHECK( alpha == Approx(1.3700514406203526e+00)); break;
        case   43: CHECK( alpha == Approx(1.2755651343706731e+00)); break;
        case   44: CHECK( alpha == Approx(1.2283219812458335e+00)); break;
        case   45: CHECK( alpha == Approx(1.2755651343706731e+00)); break;
        case   46: CHECK( alpha == Approx(1.3228082874955129e+00)); break;
        case   47: CHECK( alpha == Approx(1.5117808999948719e+00)); break;
        case   48: CHECK( alpha == Approx(1.4645377468700322e+00)); break;
        case   49: CHECK( alpha == Approx(1.4645377468700322e+00)); break;
        case   50: CHECK( alpha == Approx(1.3700514406203526e+00)); break;
        case   51: CHECK( alpha == Approx(1.3700514406203526e+00)); break;
        case   52: CHECK( alpha == Approx(1.3228082874955129e+00)); break;
        case   53: CHECK( alpha == Approx(1.3228082874955129e+00)); break;
        case   54: CHECK( alpha == Approx(1.0204521074965385e+00)); break;
        case   55: CHECK( alpha == Approx(2.5038871156165063e+00)); break;
        case   56: CHECK( alpha == Approx(2.0314555843681092e+00)); break;
        case   57: CHECK( alpha == Approx(1.8424829718687501e+00)); break;
        case   58: CHECK( alpha == Approx(1.7479966656190706e+00)); break;
        case   59: CHECK( alpha == Approx(1.7479966656190706e+00)); break;
        case   60: CHECK( alpha == Approx(1.7479966656190706e+00)); break;
        case   61: CHECK( alpha == Approx(1.7479966656190706e+00)); break;
        case   62: CHECK( alpha == Approx(1.7479966656190706e+00)); break;
        case   63: CHECK( alpha == Approx(1.7479966656190706e+00)); break;
        case   64: CHECK( alpha == Approx(1.7007535124942308e+00)); break;
        case   65: CHECK( alpha == Approx(1.6535103593693912e+00)); break;
        case   66: CHECK( alpha == Approx(1.6535103593693912e+00)); break;
        case   67: CHECK( alpha == Approx(1.6535103593693912e+00)); break;
        case   68: CHECK( alpha == Approx(1.6535103593693912e+00)); break;
        case   69: CHECK( alpha == Approx(1.6535103593693912e+00)); break;
        case   70: CHECK( alpha == Approx(1.6535103593693912e+00)); break;
        case   71: CHECK( alpha == Approx(1.6535103593693912e+00)); break;
        case   72: CHECK( alpha == Approx(1.4645377468700322e+00)); break;
        case   73: CHECK( alpha == Approx(1.3700514406203526e+00)); break;
        case   74: CHECK( alpha == Approx(1.2755651343706731e+00)); break;
        case   75: CHECK( alpha == Approx(1.2755651343706731e+00)); break;
        case   76: CHECK( alpha == Approx(1.2283219812458335e+00)); break;
        case   77: CHECK( alpha == Approx(1.2755651343706731e+00)); break;
        case   78: CHECK( alpha == Approx(1.2755651343706731e+00)); break;
        case   79: CHECK( alpha == Approx(1.2755651343706731e+00)); break;
        case   80: CHECK( alpha == Approx(1.4172945937451924e+00)); break;
        case   81: CHECK( alpha == Approx(1.7952398187439103e+00)); break;
        case   82: CHECK( alpha == Approx(1.7007535124942308e+00)); break;
        case   83: CHECK( alpha == Approx(1.5117808999948719e+00)); break;
        case   84: CHECK( alpha == Approx(1.7952398187439103e+00)); break;
        case   85: CHECK( alpha == Approx(1.1999760893709295e+00)); break;
        case   86: CHECK( alpha == Approx(1.1338356749961540e+00)); break;
        case   88: CHECK( alpha == Approx(2.0314555843681092e+00)); break;
        case   89: CHECK( alpha == Approx(1.8424829718687501e+00)); break;
        case   90: CHECK( alpha == Approx(1.7007535124942308e+00)); break;
        case   91: CHECK( alpha == Approx(1.7007535124942308e+00)); break;
        case   92: CHECK( alpha == Approx(1.6535103593693912e+00)); break;
        case   93: CHECK( alpha == Approx(1.6535103593693912e+00)); break;
        case   94: CHECK( alpha == Approx(1.6535103593693912e+00)); break;
        case   95: CHECK( alpha == Approx(1.6535103593693912e+00)); break;
        case   96: CHECK( alpha == Approx(1.8991750000000001e+00)); break;
        case   97: CHECK( alpha == Approx(1.8991750000000001e+00)); break;
        case   98: CHECK( alpha == Approx(1.8991750000000001e+00)); break;
        case   99: CHECK( alpha == Approx(1.8991750000000001e+00)); break;
        case  100: CHECK( alpha == Approx(1.8991750000000001e+00)); break;
      }
    }

  }


  SECTION("Grid Size") {
    for(auto i = 0; i < 100; ++i) {
      auto [fgr, fga] = 
        default_grid_size( AtomicNumber(i), RadialQuad::MuraKnowles,
          AtomicGridSizeDefault::FineGrid );
      REQUIRE( fgr.get() == 75 );
      REQUIRE( fga.get() == 302 );


      auto [ufgr, ufga] = 
        default_grid_size( AtomicNumber(i), RadialQuad::MuraKnowles,
          AtomicGridSizeDefault::UltraFineGrid );
      REQUIRE( ufgr.get() == 99 );
      REQUIRE( ufga.get() == 590 );
      

      auto [sfgr, sfga] = 
        default_grid_size( AtomicNumber(i), RadialQuad::MuraKnowles,
          AtomicGridSizeDefault::SuperFineGrid );
      REQUIRE( sfga.get() == 974 );
      if( i <= 2 ) REQUIRE( sfgr.get() == 175 );
      else         REQUIRE( sfgr.get() == 250 );

    }
  }

  SECTION("PySCF0 Grid Size") {
    int64_t refgr[] = { -1,
      10, 10 ,
      15, 15, 15, 15, 15, 15, 15, 15 ,
      20, 20, 20, 20, 20, 20, 20, 20,
      30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 35, 35,
      35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35,
      40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
      50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50
    };
    int64_t refga[] = { -1,
      50, 50 ,
      86, 86, 86, 86, 86, 86, 86, 86 ,
      110, 110, 110, 110, 110, 110, 110, 110,
      110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110,
      110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110,
      110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110,
      110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110
    };
    for(auto Z = 1; Z < 119; ++Z) {
      auto [gr, ga] = default_grid_size( AtomicNumber(Z), RadialQuad::MuraKnowles,
        AtomicGridSizeDefault::PySCF0 );
      REQUIRE( gr.get() == refgr[Z] );
      REQUIRE( ga.get() == refga[Z] );
    }
  }

  SECTION("PySCF1 Grid Size") {
    int64_t refgr[] = { -1,
      30, 30 ,
      40, 40, 40, 40, 40, 40, 40, 40 ,
      50, 50, 50, 50, 50, 50, 50, 50,
      60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 65, 65,
      65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65,
      70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70,
      75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75
    };
    int64_t refga[] = { -1,
      110, 110 ,
      194, 194, 194, 194, 194, 194, 194, 194 ,
      194, 194, 194, 194, 194, 194, 194, 194,
      194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194,
      194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194,
      194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194,
      194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194
    };
    for(auto Z = 1; Z < 119; ++Z) {
      auto [gr, ga] = default_grid_size( AtomicNumber(Z), RadialQuad::MuraKnowles,
        AtomicGridSizeDefault::PySCF1 );
      REQUIRE( gr.get() == refgr[Z] );
      REQUIRE( ga.get() == refga[Z] );
    }
  }

  SECTION("PySCF2 Grid Size") {
    int64_t refgr[] = { -1,
      40, 40 ,
      60, 60, 60, 60, 60, 60, 60, 60 ,
      65, 65, 65, 65, 65, 65, 65, 65,
      75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 80, 80,
      80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80,
      85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85,
      90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90
    };
    int64_t refga[] = { -1,
      194, 194 ,
      302, 302, 302, 302, 302, 302, 302, 302 ,
      302, 302, 302, 302, 302, 302, 302, 302,
      302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302,
      302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302,
      302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302,
      302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302, 302
    };
    for(auto Z = 1; Z < 119; ++Z) {
      auto [gr, ga] = default_grid_size( AtomicNumber(Z), RadialQuad::MuraKnowles,
        AtomicGridSizeDefault::PySCF2 );
      REQUIRE( gr.get() == refgr[Z] );
      REQUIRE( ga.get() == refga[Z] );
    }
  }

  SECTION("PySCF3 Grid Size") {
    int64_t refgr[] = { -1,
      50, 50 ,
      75, 75, 75, 75, 75, 75, 75, 75 ,
      80, 80, 80, 80, 80, 80, 80, 80,
      90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 95, 95,
      95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95,
      100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
      105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105
    };
    int64_t refga[] = { -1,
      302, 302 ,
      302, 302, 302, 302, 302, 302, 302, 302 ,
      434, 434, 434, 434, 434, 434, 434, 434,
      434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434,
      434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434,
      434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434,
      434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434, 434
    };
    for(auto Z = 1; Z < 119; ++Z) {
      auto [gr, ga] = default_grid_size( AtomicNumber(Z), RadialQuad::MuraKnowles,
        AtomicGridSizeDefault::PySCF3 );
      REQUIRE( gr.get() == refgr[Z] );
      REQUIRE( ga.get() == refga[Z] );
    }
  }

  SECTION("PySCF4 Grid Size") {
    int64_t refgr[] = { -1,
      60, 60 ,
      90, 90, 90, 90, 90, 90, 90, 90 ,
      95, 95, 95, 95, 95, 95, 95, 95,
      105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 110, 110,
      110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110,
      115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115,
      120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120
    };
    int64_t refga[] = { -1,
      434, 434 ,
      590, 590, 590, 590, 590, 590, 590, 590 ,
      590, 590, 590, 590, 590, 590, 590, 590,
      590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590,
      590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590,
      590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590,
      590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590, 590
    };
    for(auto Z = 1; Z < 119; ++Z) {
      auto [gr, ga] = default_grid_size( AtomicNumber(Z), RadialQuad::MuraKnowles,
        AtomicGridSizeDefault::PySCF4 );
      REQUIRE( gr.get() == refgr[Z] );
      REQUIRE( ga.get() == refga[Z] );
    }
  }

  SECTION("PySCF5 Grid Size") {
    int64_t refgr[] = { -1,
      70, 70 ,
      105, 105, 105, 105, 105, 105, 105, 105 ,
      110, 110, 110, 110, 110, 110, 110, 110,
      120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 125, 125,
      125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125,
      130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130,
      135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135
    };
    int64_t refga[] = { -1,
      590, 590 ,
      770, 770, 770, 770, 770, 770, 770, 770 ,
      770, 770, 770, 770, 770, 770, 770, 770,
      770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770,
      770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770,
      770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770,
      770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770, 770
    };
    for(auto Z = 1; Z < 119; ++Z) {
      auto [gr, ga] = default_grid_size( AtomicNumber(Z), RadialQuad::MuraKnowles,
        AtomicGridSizeDefault::PySCF5 );
      REQUIRE( gr.get() == refgr[Z] );
      REQUIRE( ga.get() == refga[Z] );
    }
  }

  SECTION("PySCF6 Grid Size") {
    int64_t refgr[] = { -1,
      80, 80 ,
      120, 120, 120, 120, 120, 120, 120, 120 ,
      125, 125, 125, 125, 125, 125, 125, 125,
      135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 140, 140,
      140, 140, 140, 140, 140, 140, 140, 140, 140, 140, 140, 140, 140, 140, 140, 140,
      145, 145, 145, 145, 145, 145, 145, 145, 145, 145, 145, 145, 145, 145, 145, 145, 145, 145, 145, 145, 145, 145, 145, 145, 145, 145, 145, 145, 145, 145, 145, 145,
      150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150
    };
    int64_t refga[] = { -1,
      770, 770 ,
      974, 974, 974, 974, 974, 974, 974, 974 ,
      974, 974, 974, 974, 974, 974, 974, 974,
      974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974,
      974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974,
      974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974,
      974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974, 974
    };
    for(auto Z = 1; Z < 119; ++Z) {
      auto [gr, ga] = default_grid_size( AtomicNumber(Z), RadialQuad::MuraKnowles,
        AtomicGridSizeDefault::PySCF6 );
      REQUIRE( gr.get() == refgr[Z] );
      REQUIRE( ga.get() == refga[Z] );
    }
  }

  SECTION("PySCF7 Grid Size") {
    int64_t refgr[] = { -1,
      90, 90 ,
      135, 135, 135, 135, 135, 135, 135, 135 ,
      140, 140, 140, 140, 140, 140, 140, 140,
      150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 155, 155,
      155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155,
      160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160,
      165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165
    };
    int64_t refga[] = { -1,
      974, 974 ,
      1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202 ,
      1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202,
      1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202,
      1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202,
      1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202,
      1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202
    };
    for(auto Z = 1; Z < 119; ++Z) {
      auto [gr, ga] = default_grid_size( AtomicNumber(Z), RadialQuad::MuraKnowles,
        AtomicGridSizeDefault::PySCF7 );
      REQUIRE( gr.get() == refgr[Z] );
      REQUIRE( ga.get() == refga[Z] );
    }
  }

  SECTION("PySCF8 Grid Size") {
    int64_t refgr[] = { -1,
      100, 100 ,
      150, 150, 150, 150, 150, 150, 150, 150 ,
      155, 155, 155, 155, 155, 155, 155, 155,
      165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 170, 170,
      170, 170, 170, 170, 170, 170, 170, 170, 170, 170, 170, 170, 170, 170, 170, 170,
      175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175,
      180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180
    };
    int64_t refga[] = { -1,
      1202, 1202 ,
      1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202 ,
      1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202,
      1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202,
      1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202,
      1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202,
      1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202, 1202
    };
    for(auto Z = 1; Z < 119; ++Z) {
      auto [gr, ga] = default_grid_size( AtomicNumber(Z), RadialQuad::MuraKnowles,
        AtomicGridSizeDefault::PySCF8 );
      REQUIRE( gr.get() == refgr[Z] );
      REQUIRE( ga.get() == refga[Z] );
    }
  }

  SECTION("PySCF9 Grid Size") {
    for(auto Z = 1; Z < 119; ++Z) {
      auto [gr, ga] = default_grid_size( AtomicNumber(Z), RadialQuad::MuraKnowles,
        AtomicGridSizeDefault::PySCF9 );
      REQUIRE( gr.get() == 200 );
      REQUIRE( ga.get() == 1454 );
    }
  }

}


TEST_CASE("Grid Specification", "[molgrid]") {

  AtomicNumber Z(6);
  auto rq  = RadialQuad::MuraKnowles;
  auto gsz = AtomicGridSizeDefault::UltraFineGrid;
  RadialSize rsz; AngularSize asz; // Clang doesn't like lambda capture of structured bindings...
  std::tie(rsz, asz) = default_grid_size(Z,rq,gsz);
  auto rscal = default_radial_scaling_factor(rq, Z);

  SECTION("Unpruned") {
    UnprunedAtomicGridSpecification gs;
    SECTION("From Sizes") {
      gs = MolGridFactory::create_default_unpruned_grid_spec(Z,rq,rsz,asz);
    }
    SECTION("From Standard") {
      gs = MolGridFactory::create_default_unpruned_grid_spec(Z,rq,gsz);
    }
    REQUIRE( gs.radial_quad  == rq    );
    REQUIRE( gs.radial_size  == rsz   );
    REQUIRE( gs.radial_scale == rscal );
    REQUIRE( gs.angular_size == asz   );
  }

  SECTION("Pruned") {

    atomic_grid_variant gs;
    std::vector<PruningRegion> ref_pruning_regions;
    UnprunedAtomicGridSpecification unp_gs =
      MolGridFactory::create_default_unpruned_grid_spec(Z,rq,gsz);
    SECTION("Unpruned") {
      gs = MolGridFactory::create_default_pruned_grid_spec(
        PruningScheme::Unpruned,Z,rq,gsz);
      ref_pruning_regions = {
        {0ul, (size_t)rsz.get(), asz}
      };
    }

    SECTION("Robust") {
      gs = MolGridFactory::create_default_pruned_grid_spec(
        PruningScheme::Robust,Z,rq,gsz);
      size_t rs = rsz.get();
      size_t r4 = rsz.get() / 4ul + 1;
      size_t r2 = rsz.get() / 2ul + 1;
      ref_pruning_regions = {
        {0ul, r4, AngularSize(26)},
        {r4,  r2, AngularSize(434)},
        {r2,  rs, AngularSize(590)}
      };
    }

    SECTION("Treutler") {
      gs = MolGridFactory::create_default_pruned_grid_spec(
        PruningScheme::Treutler,Z,rq,gsz);
      size_t rs = rsz.get();
      size_t r3 = rsz.get() / 3ul + 1;
      size_t r2 = rsz.get() / 2ul + 1;
      ref_pruning_regions = {
        {0ul, r3, AngularSize(26)},
        {r3,  r2, AngularSize(50)},
        {r2,  rs, AngularSize(590)}
      };
    }

    SECTION("PySCF Treutler") {
      gs = MolGridFactory::create_default_pruned_grid_spec(
        PruningScheme::PySCF_Treutler,Z,rq,gsz);
      size_t rs = rsz.get();
      size_t r3 = rsz.get() / 3ul;
      size_t r2 = rsz.get() / 2ul;
      ref_pruning_regions = {
        {0ul, r3, AngularSize(14)},
        {r3,  r2, AngularSize(50)},
        {r2,  rs, AngularSize(590)}
      };
    }

    std::visit([=](auto& g) {
      REQUIRE( g.radial_quad  == rq    );
      REQUIRE( g.radial_size  == rsz   );
      REQUIRE( g.radial_scale == rscal );
    }, gs);

    // Pruned Check
    if( const auto* pru = std::get_if<PrunedAtomicGridSpecification>(&gs) ) {
      REQUIRE( pru->pruning_regions == ref_pruning_regions );
    } 

  }

  SECTION("PySCF Slater-Bragg Radii") {

    const std::vector<std::pair<int,double>> pyscf_slater_reference = {
      {  1, 6.61404143597771554e-01},
      {  2, 2.64561657439108622e+00},
      {  6, 1.32280828719554311e+00},
      { 10, 2.83458918684759276e+00},
      { 18, 3.40150702421711149e+00},
      { 36, 3.59047963667361714e+00},
      { 54, 3.96842486158663021e+00},
      { 86, 3.96842486158663021e+00},
      {103, 3.30702071798885822e+00},
      {130, 3.30702071798885822e+00},
    };

    for( const auto& [Z, radius] : pyscf_slater_reference ) {
      CHECK( pyscf_slater_radius_64(AtomicNumber(Z)) == Approx(radius) );
    }

  }

}


TEST_CASE("PySCF Pruning Specifications", "[molgrid]") {

  struct PySCFPruningCase {
    PruningScheme scheme;
    AtomicNumber Z;
    AtomicGridSizeDefault grid_size;
    std::vector<PruningRegion> reference_regions;
  };

  const std::vector<PySCFPruningCase> test_cases = {
    {PruningScheme::PySCF_Treutler, AtomicNumber(1), AtomicGridSizeDefault::PySCF0,
      {{0ul, 3ul, AngularSize(14)}, {3ul, 10ul, AngularSize(50)}}},
    {PruningScheme::PySCF_Treutler, AtomicNumber(11), AtomicGridSizeDefault::PySCF3,
      {{0ul, 26ul, AngularSize(14)}, {26ul, 40ul, AngularSize(50)},
       {40ul, 80ul, AngularSize(434)}}},
    {PruningScheme::PySCF_Treutler, AtomicNumber(8), AtomicGridSizeDefault::PySCF2,
      {{0ul, 20ul, AngularSize(14)}, {20ul, 30ul, AngularSize(50)},
       {30ul, 60ul, AngularSize(302)}}},
    {PruningScheme::PySCF_Treutler, AtomicNumber(36), AtomicGridSizeDefault::PySCF4,
      {{0ul, 35ul, AngularSize(14)}, {35ul, 52ul, AngularSize(50)},
       {52ul, 105ul, AngularSize(590)}}},

    {PruningScheme::PySCF_SG1, AtomicNumber(1), AtomicGridSizeDefault::PySCF0,
      {{0ul, 3ul, AngularSize(6)}, {3ul, 4ul, AngularSize(38)},
       {4ul, 5ul, AngularSize(86)}, {5ul, 9ul, AngularSize(194)},
       {9ul, 10ul, AngularSize(86)}}},
    {PruningScheme::PySCF_SG1, AtomicNumber(2), AtomicGridSizeDefault::PySCF4,
      {{0ul, 17ul, AngularSize(6)}, {17ul, 21ul, AngularSize(38)},
       {21ul, 26ul, AngularSize(86)}, {26ul, 42ul, AngularSize(194)},
       {42ul, 60ul, AngularSize(86)}}},
    {PruningScheme::PySCF_SG1, AtomicNumber(6), AtomicGridSizeDefault::PySCF3,
      {{0ul, 22ul, AngularSize(6)}, {22ul, 31ul, AngularSize(38)},
       {31ul, 38ul, AngularSize(86)}, {38ul, 57ul, AngularSize(194)},
       {57ul, 75ul, AngularSize(86)}}},
    {PruningScheme::PySCF_SG1, AtomicNumber(8), AtomicGridSizeDefault::PySCF2,
      {{0ul, 17ul, AngularSize(6)}, {17ul, 24ul, AngularSize(38)},
       {24ul, 29ul, AngularSize(86)}, {29ul, 44ul, AngularSize(194)},
       {44ul, 60ul, AngularSize(86)}}},
    {PruningScheme::PySCF_SG1, AtomicNumber(11), AtomicGridSizeDefault::PySCF1,
      {{0ul, 17ul, AngularSize(6)}, {17ul, 26ul, AngularSize(38)},
       {26ul, 33ul, AngularSize(86)}, {33ul, 45ul, AngularSize(194)},
       {45ul, 50ul, AngularSize(86)}}},
    {PruningScheme::PySCF_SG1, AtomicNumber(18), AtomicGridSizeDefault::PySCF4,
      {{0ul, 25ul, AngularSize(6)}, {25ul, 39ul, AngularSize(38)},
       {39ul, 48ul, AngularSize(86)}, {48ul, 69ul, AngularSize(194)},
       {69ul, 95ul, AngularSize(86)}}},

    {PruningScheme::PySCF_NWChem, AtomicNumber(1), AtomicGridSizeDefault::PySCF0,
      {{0ul, 3ul, AngularSize(50)}, {3ul, 8ul, AngularSize(74)},
       {8ul, 10ul, AngularSize(50)}}},
    {PruningScheme::PySCF_NWChem, AtomicNumber(2), AtomicGridSizeDefault::PySCF4,
      {{0ul, 27ul, AngularSize(50)}, {27ul, 34ul, AngularSize(86)},
       {34ul, 42ul, AngularSize(350)}, {42ul, 59ul, AngularSize(434)},
       {59ul, 60ul, AngularSize(350)}}},
    {PruningScheme::PySCF_NWChem, AtomicNumber(6), AtomicGridSizeDefault::PySCF3,
      {{0ul, 22ul, AngularSize(50)}, {22ul, 32ul, AngularSize(86)},
       {32ul, 38ul, AngularSize(266)}, {38ul, 58ul, AngularSize(302)},
       {58ul, 75ul, AngularSize(266)}}},
    {PruningScheme::PySCF_NWChem, AtomicNumber(8), AtomicGridSizeDefault::PySCF2,
      {{0ul, 18ul, AngularSize(50)}, {18ul, 26ul, AngularSize(86)},
       {26ul, 31ul, AngularSize(266)}, {31ul, 47ul, AngularSize(302)},
       {47ul, 60ul, AngularSize(266)}}},
    {PruningScheme::PySCF_NWChem, AtomicNumber(10), AtomicGridSizeDefault::PySCF3,
      {{0ul, 31ul, AngularSize(50)}, {31ul, 43ul, AngularSize(86)},
       {43ul, 52ul, AngularSize(266)}, {52ul, 72ul, AngularSize(302)},
       {72ul, 75ul, AngularSize(266)}}},
    {PruningScheme::PySCF_NWChem, AtomicNumber(36), AtomicGridSizeDefault::PySCF4,
      {{0ul, 39ul, AngularSize(50)}, {39ul, 61ul, AngularSize(86)},
       {61ul, 75ul, AngularSize(434)}, {75ul, 99ul, AngularSize(590)},
       {99ul, 105ul, AngularSize(434)}}},

    {PruningScheme::PySCF_SGX, AtomicNumber(1), AtomicGridSizeDefault::PySCF5,
      {{0ul, 21ul, AngularSize(50)}, {21ul, 26ul, AngularSize(302)},
       {26ul, 33ul, AngularSize(434)}, {33ul, 52ul, AngularSize(590)},
       {52ul, 70ul, AngularSize(434)}}},
    {PruningScheme::PySCF_SGX, AtomicNumber(2), AtomicGridSizeDefault::PySCF0,
      {{0ul, 4ul, AngularSize(14)}, {4ul, 6ul, AngularSize(26)},
       {6ul, 10ul, AngularSize(50)}}},
    {PruningScheme::PySCF_SGX, AtomicNumber(6), AtomicGridSizeDefault::PySCF4,
      {{0ul, 27ul, AngularSize(50)}, {27ul, 38ul, AngularSize(302)},
       {38ul, 46ul, AngularSize(434)}, {46ul, 70ul, AngularSize(590)},
       {70ul, 90ul, AngularSize(434)}}},
    {PruningScheme::PySCF_SGX, AtomicNumber(8), AtomicGridSizeDefault::PySCF2,
      {{0ul, 18ul, AngularSize(26)}, {18ul, 26ul, AngularSize(110)},
       {26ul, 31ul, AngularSize(194)}, {31ul, 47ul, AngularSize(302)},
       {47ul, 60ul, AngularSize(194)}}},
    {PruningScheme::PySCF_SGX, AtomicNumber(10), AtomicGridSizeDefault::PySCF3,
      {{0ul, 31ul, AngularSize(26)}, {31ul, 43ul, AngularSize(110)},
       {43ul, 52ul, AngularSize(194)}, {52ul, 72ul, AngularSize(302)},
       {72ul, 75ul, AngularSize(194)}}},
    {PruningScheme::PySCF_SGX, AtomicNumber(18), AtomicGridSizeDefault::PySCF4,
      {{0ul, 34ul, AngularSize(50)}, {34ul, 52ul, AngularSize(302)},
       {52ul, 65ul, AngularSize(434)}, {65ul, 87ul, AngularSize(590)},
       {87ul, 95ul, AngularSize(434)}}}
  };

  for( const auto& test_case : test_cases ) {
    auto gs = MolGridFactory::create_default_pruned_grid_spec(
      test_case.scheme, test_case.Z, RadialQuad::TreutlerAhlrichs,
      test_case.grid_size
    );

    const auto* pruned = std::get_if<PrunedAtomicGridSpecification>(&gs);
    REQUIRE( pruned );
    CHECK( pruned->pruning_regions == test_case.reference_regions );
  }

}

TEST_CASE("PySCF Pruning Known Differences", "[molgrid]") {

  SECTION("PySCF Treutler radial grids extend past GauXC TA defaults") {
    REQUIRE_THROWS( MolGridFactory::create_default_pruned_grid_spec(
      PruningScheme::PySCF_NWChem, AtomicNumber(54),
      RadialQuad::TreutlerAhlrichs, AtomicGridSizeDefault::PySCF4
    ) );
  }

}

TEST_CASE("PySCF Pruning Exceptions", "[molgrid]") {

  const UnprunedAtomicGridSpecification unp{
    RadialQuad::TreutlerAhlrichs, RadialSize(10), RadialScale(1.0),
    AngularSize(302)
  };

  SECTION("Atom-aware PySCF schemes require an atomic number") {
    REQUIRE_THROWS_AS( create_pruned_spec(PruningScheme::PySCF_SG1, unp),
                       generic_gauxc_exception );
    REQUIRE_THROWS_AS( create_pruned_spec(PruningScheme::PySCF_NWChem, unp),
                       generic_gauxc_exception );
    REQUIRE_THROWS_AS( create_pruned_spec(PruningScheme::PySCF_SGX, unp),
                       generic_gauxc_exception );
  }

  SECTION("SG1 radii are limited to H-Ar; PySCF raises IndexError for Z=19") {
    REQUIRE_THROWS_AS(
      create_pruned_spec(PruningScheme::PySCF_SG1, AtomicNumber(19), unp),
      generic_gauxc_exception
    );
  }

  SECTION("PySCF accepts Z=0 SG1 via its dummy radius; GauXC rejects it") {
    REQUIRE_THROWS_AS(
      create_pruned_spec(PruningScheme::PySCF_SG1, AtomicNumber(0), unp),
      generic_gauxc_exception
    );
  }

  SECTION("Bragg radii stop at Z=130; PySCF raises IndexError for Z=131") {
    REQUIRE_THROWS_AS(
      create_pruned_spec(PruningScheme::PySCF_NWChem, AtomicNumber(131), unp),
      generic_gauxc_exception
    );
    REQUIRE_THROWS_AS(
      create_pruned_spec(PruningScheme::PySCF_SGX, AtomicNumber(131), unp),
      generic_gauxc_exception
    );
  }

  SECTION("NWChem requires a supported Lebedev angular size; PySCF raises IndexError") {
    auto invalid_angular = unp;
    invalid_angular.angular_size = AngularSize(52);
    REQUIRE_THROWS_AS(
      create_pruned_spec(PruningScheme::PySCF_NWChem, AtomicNumber(8), invalid_angular),
      generic_gauxc_exception
    );
  }

  SECTION("SGX is bounded by PySCF's mapping table; PySCF raises IndexError") {
    auto too_large_sgx = unp;
    too_large_sgx.angular_size = AngularSize(6000);
    REQUIRE_THROWS_AS(
      create_pruned_spec(PruningScheme::PySCF_SGX, AtomicNumber(8), too_large_sgx),
      generic_gauxc_exception
    );
  }

}



#if 0

TEST_CASE("MolGrid", "[molgrid]") {

  // Water
  Molecule mol;
  mol.emplace_back(AtomicNumber(1), 0., 1.579252144093028,  2.174611055780858);
  mol.emplace_back(AtomicNumber(8), 0., 0.000000000000000,  0.000000000000000);
  mol.emplace_back(AtomicNumber(1), 0., 1.579252144093028, -2.174611055780858);

  auto mk_scaling = 
    default_radial_scaling_factors( RadialQuad::MuraKnowles, AtomicNumber(8) );
  auto mhl_scaling = 
    default_radial_scaling_factors( RadialQuad::MurrayHandyLaming, AtomicNumber(8) );

  auto [fg_radial_size, fg_angular_size] =
    default_grid_size( AtomicGridSizeDefault::FineGrid, AtomicNumber(8) );
  auto [ufg_radial_size, ufg_angular_size] =
    default_grid_size( AtomicGridSizeDefault::UltraFineGrid, AtomicNumber(8) );
  auto [sfg_radial_size, sfg_angular_size] =
    default_grid_size( AtomicGridSizeDefault::SuperFineGrid, AtomicNumber(8) );

  SECTION("MK") {

    SECTION("Explicit Construction") {

      MolGrid mg( RadialQuad::MuraKnowles, fg_size, mk_scaling, mol );
      
      CHECK( mg.natoms_uniq() == 2 );

      for( const auto& atom : mol ) {
        auto alpha = mk_scaling[atom.Z];
        auto gsz   = fg_size[atom.Z];

        Grid atom_grid( RadialQuad::MuraKnowles, gsz, alpha );

        //CHECK( mg.get_rscal_factor(atom.Z) == alpha );
        CHECK( mg.get_grid_size( atom. Z ) == gsz   );

        CHECK( atom_grid.batcher().quadrature().points() ==
               mg.get_grid(atom.Z).batcher().quadrature().points() );
        CHECK( atom_grid.batcher().quadrature().weights() ==
               mg.get_grid(atom.Z).batcher().quadrature().weights() );

      }

    }

    SECTION("Default Scaling Factors") {

      MolGrid mg( RadialQuad::MuraKnowles, fg_size, mol );
      
      CHECK( mg.natoms_uniq() == 2 );

      for( const auto& atom : mol ) {
        auto alpha = mk_scaling[atom.Z];
        auto gsz   = fg_size[atom.Z];

        Grid atom_grid( RadialQuad::MuraKnowles, gsz, alpha );

        //CHECK( mg.get_rscal_factor(atom.Z) == alpha );
        CHECK( mg.get_grid_size( atom. Z ) == gsz   );

        CHECK( atom_grid.batcher().quadrature().points() ==
               mg.get_grid(atom.Z).batcher().quadrature().points() );
        CHECK( atom_grid.batcher().quadrature().weights() ==
               mg.get_grid(atom.Z).batcher().quadrature().weights() );

      }

    }

    SECTION("Named Default Grid Size") {

      MolGrid mg( RadialQuad::MuraKnowles, AtomicGridSizeDefault::UltraFineGrid, 
        mk_scaling, mol );
      
      CHECK( mg.natoms_uniq() == 2 );

      for( const auto& atom : mol ) {
        auto alpha = mk_scaling[atom.Z];
        auto gsz   = ufg_size[atom.Z];

        Grid atom_grid( RadialQuad::MuraKnowles, gsz, alpha );

        //CHECK( mg.get_rscal_factor(atom.Z) == alpha );
        CHECK( mg.get_grid_size( atom. Z ) == gsz   );

        CHECK( atom_grid.batcher().quadrature().points() ==
               mg.get_grid(atom.Z).batcher().quadrature().points() );
        CHECK( atom_grid.batcher().quadrature().weights() ==
               mg.get_grid(atom.Z).batcher().quadrature().weights() );

      }

    }

    SECTION("Named Default Grid Size + Default Scaling") {

      MolGrid mg( RadialQuad::MuraKnowles, AtomicGridSizeDefault::SuperFineGrid, 
        mol );
      
      CHECK( mg.natoms_uniq() == 2 );

      for( const auto& atom : mol ) {
        auto alpha = mk_scaling[atom.Z];
        auto gsz   = sfg_size[atom.Z];

        Grid atom_grid( RadialQuad::MuraKnowles, gsz, alpha );

        //CHECK( mg.get_rscal_factor(atom.Z) == alpha );
        CHECK( mg.get_grid_size( atom. Z ) == gsz   );

        CHECK( atom_grid.batcher().quadrature().points() ==
               mg.get_grid(atom.Z).batcher().quadrature().points() );
        CHECK( atom_grid.batcher().quadrature().weights() ==
               mg.get_grid(atom.Z).batcher().quadrature().weights() );

      }

    }

  }

  SECTION("MHL") {

    MolGrid mg( RadialQuad::MurrayHandyLaming, fg_size, mhl_scaling, mol );
    
    CHECK( mg.natoms_uniq() == 2 );

    for( const auto& atom : mol ) {
      auto alpha = mhl_scaling[atom.Z];
      auto gsz   = fg_size[atom.Z];

      Grid atom_grid( RadialQuad::MurrayHandyLaming, gsz, alpha );

      //CHECK( mg.get_rscal_factor(atom.Z) == alpha );
      CHECK( mg.get_grid_size( atom. Z ) == gsz   );

      CHECK( atom_grid.batcher().quadrature().points() ==
             mg.get_grid(atom.Z).batcher().quadrature().points() );
      CHECK( atom_grid.batcher().quadrature().weights() ==
             mg.get_grid(atom.Z).batcher().quadrature().weights() );

    }

  }

#if 0
  SECTION("Default") {
    MolGrid mg( AtomicGridSizeDefault::FineGrid, mol );
    for( const auto& atom: mol )
      CHECK( mg.get_radial_quad(atom.Z) == RadialQuad::MuraKnowles );
  }
#endif

}
#endif



