#include <gauxc/molgrid/defaults.hpp>
#include <gauxc/exceptions.hpp>
#include <integratorxx/quadratures/lebedev_laikov.hpp>

namespace GauXC {

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
  const double fac = (Z==1) ? 0.5 : 1.0;
  return RadialScale( default_atomic_radius(_Z) * fac );
}

RadialScale default_radial_scaling_factor(RadialQuad rq, AtomicNumber Z) {
  if( rq == RadialQuad::MuraKnowles ) 
    return default_mk_radial_scaling_factor(Z);
  else if( rq == RadialQuad::TreutlerAldrichs )
    return default_ta_radial_scaling_factor(Z);
  else // MHL
    return default_mhl_radial_scaling_factor(Z);
}

PrunedAtomicGridSpecification robust_psi4_pruning_scheme(
  UnprunedAtomicGridSpecification unp ) {

  // Look up order
  // XXX: THIS ONLY WORKS FOR LEBEDEV
  using namespace IntegratorXX::detail::lebedev;
  const auto asz = unp.angular_size.get();
  const auto base_order = algebraic_order_by_npts(asz);
  if( base_order < 0 ) GAUXC_GENERIC_EXCEPTION("Invalid Base Grid");

  const auto med_order = 
    next_algebraic_order(base_order > 6 ? base_order-6 : base_order);
  const auto low_order = 7;

  AngularSize med_sz(npts_by_algebraic_order(med_order));
  AngularSize low_sz(npts_by_algebraic_order(low_order));

  // Create Pruning Regions
  const size_t rsz = unp.radial_size.get();
  const size_t r_div_4 = rsz / 4ul;
  std::vector<PruningRegion> pruning_regions = {
    {  0ul,       r_div_4,   low_sz},
    {  r_div_4, 2ul*r_div_4, med_sz},
    {2ul*r_div_4,       rsz, unp.angular_size}
  };

  return PrunedAtomicGridSpecification{
    unp.radial_quad, unp.radial_size, unp.radial_scale, pruning_regions
  };
  
}




#if 0

atomic_scal_factor_map get_default_scaling_factors( RadialQuad rq, 
  AtomicNumber maxZ ) {


  atomic_scal_factor_map smap;

  if( rq == RadialQuad::MuraKnowles ) {

    for( int Z = 1; Z <= maxZ.get(); ++Z ) {

      switch (Z) {

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
          smap.insert_or_assign( AtomicNumber(Z), RadialScale(7.0) );
          break;
        default:
          smap.insert_or_assign( AtomicNumber(Z), RadialScale(5.0) );
          break;

      }

    }

  } else if( rq == RadialQuad::TreutlerAldrichs ) {

    smap = {
      { AtomicNumber(1),  RadialScale(0.8) }, // H
      { AtomicNumber(2),  RadialScale(0.9) }, // He
      { AtomicNumber(3),  RadialScale(1.8) }, // Li
      { AtomicNumber(4),  RadialScale(1.4) }, // Be
      { AtomicNumber(5),  RadialScale(1.3) }, // B
      { AtomicNumber(6),  RadialScale(1.1) }, // C
      { AtomicNumber(7),  RadialScale(0.9) }, // N
      { AtomicNumber(8),  RadialScale(0.9) }, // O
      { AtomicNumber(9),  RadialScale(0.9) }, // F
      { AtomicNumber(10), RadialScale(0.9) }, // Ne
      { AtomicNumber(11), RadialScale(1.4) }, // Na
      { AtomicNumber(12), RadialScale(1.3) }, // Mg
      { AtomicNumber(13), RadialScale(1.3) }, // Al
      { AtomicNumber(14), RadialScale(1.2) }, // Si
      { AtomicNumber(15), RadialScale(1.1) }, // P
      { AtomicNumber(16), RadialScale(1.0) }, // S
      { AtomicNumber(17), RadialScale(1.0) }, // Cl
      { AtomicNumber(18), RadialScale(1.0) }, // Ar
      { AtomicNumber(19), RadialScale(1.5) }, // K
      { AtomicNumber(20), RadialScale(1.4) }, // Ca
      { AtomicNumber(21), RadialScale(1.3) }, // Sc
      { AtomicNumber(22), RadialScale(1.2) }, // Ti
      { AtomicNumber(23), RadialScale(1.2) }, // V
      { AtomicNumber(24), RadialScale(1.2) }, // Cr
      { AtomicNumber(25), RadialScale(1.2) }, // Mn
      { AtomicNumber(26), RadialScale(1.2) }, // Fe
      { AtomicNumber(27), RadialScale(1.2) }, // Co
      { AtomicNumber(28), RadialScale(1.1) }, // Ni
      { AtomicNumber(29), RadialScale(1.1) }, // Cu
      { AtomicNumber(30), RadialScale(1.1) }, // Zn
      { AtomicNumber(31), RadialScale(1.1) }, // Ga
      { AtomicNumber(32), RadialScale(1.0) }, // Ge
      { AtomicNumber(33), RadialScale(0.9) }, // As
      { AtomicNumber(34), RadialScale(0.9) }, // Se
      { AtomicNumber(35), RadialScale(0.9) }, // Br
      { AtomicNumber(36), RadialScale(0.9) }  // Kr
    };
    

  } else { // Slater Radii + missing terms through

    auto slater_64 = slater_radii_64();
    auto clementi_67 = clementi_radii_67();

    smap = slater_64;
    smap.merge( clementi_67 );

    auto maxRefZ = std::max_element( smap.begin(), smap.end(),
      [](const auto& a, const auto &b) {
        return a.first.get() < b.first.get();
      })->first;

    // Fill in remaining with 2.01 Angstroms -> 3.79825 Bohr
    if( maxRefZ.get() < maxZ.get() )
    for( int64_t i = maxRefZ.get() + 1; i <= maxZ.get(); ++i )
      smap.insert_or_assign( AtomicNumber(i), RadialScale(3.79835) );

    for( auto& [Z, R] : smap )
      if( Z != 1ll ) R = RadialScale( R.get() / 2. );

  }


  return smap;
}




atomic_grid_size_map get_default_grid_sizes( AtomicGridSizeDefault G, 
  AtomicNumber maxZ ) {

  atomic_grid_size_map map;

  if( G == AtomicGridSizeDefault::FineGrid ) {

    auto rsz = RadialSize(75);
    auto asz = AngularSize(302);
    for( int64_t i = 1; i <= maxZ.get(); ++i )
      map.insert_or_assign( AtomicNumber(i), std::pair{ rsz, asz } );

  } else if( G == AtomicGridSizeDefault::UltraFineGrid ) {

    auto rsz = RadialSize(99);
    auto asz = AngularSize(590);
    for( int64_t i = 1; i <= maxZ.get(); ++i )
      map.insert_or_assign( AtomicNumber(i), std::pair{ rsz, asz } );

  } else if( G == AtomicGridSizeDefault::SuperFineGrid ) {

    auto rsz_sm = RadialSize(175);
    auto rsz_bg = RadialSize(250);
    auto asz = AngularSize(974);

    for( int64_t i = 1; i <= maxZ.get(); ++i ) {
      if( i <= 2 )
        map.insert_or_assign( AtomicNumber(i), std::pair{ rsz_sm, asz } );
      else
        map.insert_or_assign( AtomicNumber(i), std::pair{ rsz_bg, asz } );
    }

  }

  return map;
}


#endif
}
