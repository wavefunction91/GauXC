#pragma once

#include <mpi.h>
#include <memory>

#include "types.hpp"


namespace GauXC {

class XCIntegratorBase {

protected:

  MPI_Comm         comm_;
  functional_type  func_;
  BasisSet         basis_;
  Molecule         mol_;
  
  
public:

  XCIntegratorBase() = delete;

  XCIntegratorBase( MPI_Comm comm, const functional_type& func,
    const BasisSet& basis, const Molecule& mol );

  XCIntegratorBase( const XCIntegratorBase& ) = default;
  XCIntegratorBase( XCIntegratorBase&& )      = default;


  void genseg_quad( RadialQuad rquad, int64_t n_rad, int64_t n_ang );


};

template <typename MatrixType>
class XCIntegrator : public XCIntegratorBase {

public:

  using value_type = typename MatrixType::value_type;  



};

}
