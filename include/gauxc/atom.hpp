#pragma once

#include "types.hpp"

namespace GauXC {

using AtomicNumber = detail::NamedType< int64_t, struct AtomicNumberType >;

struct Atom {

  AtomicNumber Z;

  double x;
  double y;
  double z;

  Atom( AtomicNumber _Z, double _x, double _y, double _z ) :
    Z(_Z), x(_x), y(_y), z(_z) { }

};

inline bool operator==( const Atom& a1, const Atom& a2 ) {
  return a1.Z == a2.Z and a1.x == a2.x and a1.y == a2.y and a1.z == a2.z; 
}

}
