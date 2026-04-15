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
#include <gauxc/molgrid/defaults.hpp>
#include <gauxc/exceptions.hpp>
#include <integratorxx/quadratures/s2/lebedev_laikov.hpp>

namespace GauXC {

namespace detail {

inline RadialSize get_pyscf_radial_size( AtomicNumber Z, int level ) {
  if (level < 0 || level > 8) 
    GAUXC_GENERIC_EXCEPTION("Invalid PySCF grid level: " + std::to_string(level) + ". Valid levels are 0-8.");

  if ( Z.get() <= 2 ) {
    if (level == 0) return RadialSize( 10 );
    return RadialSize( 20 + 10*level );  // ..., 30, 40, 50, ...
  } else if ( Z.get() <= 10 ) {
    if (level == 0) return RadialSize( 15 );
    if (level == 1) return RadialSize( 40 );
    return RadialSize( 30 + 15 * level );  // ..., 60, 75, 90, ...
  } else if ( Z.get() <= 18 ) {
    if (level == 0) return RadialSize( 20 );
    return RadialSize( 35 + 15 * level );  // ..., 50, 65, 80, ...
  } else if ( Z.get() <= 36 ) {
    if (level == 0) return RadialSize( 30 );
    return RadialSize( 45 + 15 * level );  // ..., 60, 75, 90, ...
  } else if ( Z.get() <= 54 ) {
    if (level == 0) return RadialSize( 35 );
    return RadialSize( 50 + 15 * level );  // ..., 65, 80, 95, ...
  } else if ( Z.get() <= 86 ) {
    if (level == 0) return RadialSize( 40 );
    return RadialSize( 55 + 15 * level );  // ..., 70, 85, 100, ...
  } else if ( Z.get() <= 118 ) {
    if (level == 0) return RadialSize( 50 );
    return RadialSize( 60 + 15 * level );  // ..., 75, 90, 105, ...
  } else {
    GAUXC_GENERIC_EXCEPTION("Z > 118 Not Supported for PySCF Grid Defaults");
  }
}
}

RadialScale default_mk_radial_scaling_factor( AtomicNumber _Z ) {
  auto Z = _Z.get();
  switch(Z) {
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
      return RadialScale(7.0);
    default:
      return RadialScale(5.0);
  }
}

RadialScale default_ta_radial_scaling_factor( AtomicNumber _Z ) {
  auto Z = _Z.get();
  switch(Z) {
    case 1: return  RadialScale(0.8); // H
    case 2: return  RadialScale(0.9); // He
    case 3: return  RadialScale(1.8); // Li
    case 4: return  RadialScale(1.4); // Be
    case 5: return  RadialScale(1.3); // B
    case 6: return  RadialScale(1.1); // C
    case 7: return  RadialScale(0.9); // N
    case 8: return  RadialScale(0.9); // O
    case 9: return  RadialScale(0.9); // F
    case 10: return RadialScale(0.9); // Ne
    case 11: return RadialScale(1.4); // Na
    case 12: return RadialScale(1.3); // Mg
    case 13: return RadialScale(1.3); // Al
    case 14: return RadialScale(1.2); // Si
    case 15: return RadialScale(1.1); // P
    case 16: return RadialScale(1.0); // S
    case 17: return RadialScale(1.0); // Cl
    case 18: return RadialScale(1.0); // Ar
    case 19: return RadialScale(1.5); // K
    case 20: return RadialScale(1.4); // Ca
    case 21: return RadialScale(1.3); // Sc
    case 22: return RadialScale(1.2); // Ti
    case 23: return RadialScale(1.2); // V
    case 24: return RadialScale(1.2); // Cr
    case 25: return RadialScale(1.2); // Mn
    case 26: return RadialScale(1.2); // Fe
    case 27: return RadialScale(1.2); // Co
    case 28: return RadialScale(1.1); // Ni
    case 29: return RadialScale(1.1); // Cu
    case 30: return RadialScale(1.1); // Zn
    case 31: return RadialScale(1.1); // Ga
    case 32: return RadialScale(1.0); // Ge
    case 33: return RadialScale(0.9); // As
    case 34: return RadialScale(0.9); // Se
    case 35: return RadialScale(0.9); // Br
    case 36: return RadialScale(0.9); // Kr
    default:
      GAUXC_GENERIC_EXCEPTION("Z > 36 Not Supported for TA Quadrature");
      abort();
  }
}

RadialScale default_mhl_radial_scaling_factor( AtomicNumber _Z ) {
  auto Z = _Z.get(); 
  const double fac = (Z!=1) ? 0.5 : 1.0;
  return RadialScale( default_atomic_radius(_Z) * fac );
}

RadialScale default_bk_radial_scaling_factor( AtomicNumber _Z ) {
  auto Z = _Z.get(); 
  const double fac = (Z!=1) ? 0.5 : 1.0;
  return RadialScale( default_atomic_radius(_Z) * fac );
}

RadialScale default_radial_scaling_factor(RadialQuad rq, AtomicNumber Z) {
  if( rq == RadialQuad::MuraKnowles ) 
    return default_mk_radial_scaling_factor(Z);
  else if( rq == RadialQuad::TreutlerAhlrichs )
    return default_ta_radial_scaling_factor(Z);
  else if( rq == RadialQuad::Becke )
    return default_bk_radial_scaling_factor(Z);
  else // MHL
    return default_mhl_radial_scaling_factor(Z);
}




std::tuple<RadialSize,AngularSize> 
  default_grid_size(AtomicNumber Z, RadialQuad /*rq*/, AtomicGridSizeDefault s) {

  switch(s) {
    case AtomicGridSizeDefault::GM3:
      return std::make_tuple( RadialSize(35), AngularSize(110) );

    case AtomicGridSizeDefault::GM5:
      return std::make_tuple( RadialSize(50), AngularSize(302) );

    case AtomicGridSizeDefault::PySCF0:
      return std::make_tuple( detail::get_pyscf_radial_size(Z, 0), Z.get() <= 2 ? AngularSize(50) : (Z.get() <= 10 ? AngularSize(86) : AngularSize(110)));

    case AtomicGridSizeDefault::PySCF1:
      return std::make_tuple( detail::get_pyscf_radial_size(Z, 1), Z.get() <= 2 ? AngularSize(110) : AngularSize(194));

    case AtomicGridSizeDefault::PySCF2:
      return std::make_tuple( detail::get_pyscf_radial_size(Z, 2), Z.get() <= 2 ? AngularSize(194) : AngularSize(302));

    case AtomicGridSizeDefault::PySCF3:
      return std::make_tuple( detail::get_pyscf_radial_size(Z, 3), Z.get() <= 10 ? AngularSize(302) : AngularSize(434));

    case AtomicGridSizeDefault::PySCF4:
      return std::make_tuple( detail::get_pyscf_radial_size(Z, 4), Z.get() <= 2 ? AngularSize(434) : AngularSize(590));

    case AtomicGridSizeDefault::PySCF5:
      return std::make_tuple( detail::get_pyscf_radial_size(Z, 5), Z.get() <= 2 ? AngularSize(590) : AngularSize(770));

    case AtomicGridSizeDefault::PySCF6:
      return std::make_tuple( detail::get_pyscf_radial_size(Z, 6), Z.get() <= 2 ? AngularSize(770) : AngularSize(974));

    case AtomicGridSizeDefault::PySCF7:
      return std::make_tuple( detail::get_pyscf_radial_size(Z, 7), Z.get() <= 2 ? AngularSize(974) : AngularSize(1202));

    case AtomicGridSizeDefault::PySCF8:
      return std::make_tuple( detail::get_pyscf_radial_size(Z, 8), AngularSize(1202));

    case AtomicGridSizeDefault::PySCF9:
      return std::make_tuple( RadialSize(200), AngularSize(1454) );

    case AtomicGridSizeDefault::FineGrid:
      return std::make_tuple( RadialSize(75), AngularSize(302) );

    case AtomicGridSizeDefault::UltraFineGrid:
      return std::make_tuple( RadialSize(99), AngularSize(590) );

    case AtomicGridSizeDefault::SuperFineGrid:
      if( Z.get() <= 2 ) {
        return std::make_tuple( RadialSize(175), AngularSize(974) );
      } else {
        return std::make_tuple( RadialSize(250), AngularSize(974) );
      }

    default:
      GAUXC_GENERIC_EXCEPTION("Not A Recognized Standard Grid");
      abort();
  }

}



UnprunedAtomicGridSpecification MolGridFactory::create_default_unpruned_grid_spec(
  AtomicNumber Z, RadialQuad rq, RadialSize rsz, AngularSize asz
) {
  return UnprunedAtomicGridSpecification{
    rq, rsz, default_radial_scaling_factor(rq,Z), asz
  };
}

UnprunedAtomicGridSpecification MolGridFactory::create_default_unpruned_grid_spec(
  AtomicNumber Z, RadialQuad rq, AtomicGridSizeDefault standard_grid
) {
  auto [rsz, asz] = default_grid_size(Z, rq, standard_grid);
  return create_default_unpruned_grid_spec(Z, rq, rsz, asz);
}





}
