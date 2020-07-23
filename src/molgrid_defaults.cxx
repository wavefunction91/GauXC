#include <gauxc/molgrid/defaults.hpp>

namespace GauXC {


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



}
