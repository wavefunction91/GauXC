#pragma once

#include <mpi.h>
#include <memory>

#include <gauxc/types.hpp>
#include <gauxc/basisset.hpp>
#include <gauxc/molecule.hpp>
#include <gauxc/molgrid.hpp>

namespace GauXC {

namespace detail {
  class XCIntegratorImpl;
}


template <typename MatrixType>
class XCIntegrator : public XCIntegratorBase {

public:

  using matrix_type = MatrixType;
  using value_type  = typename matrix_type::value_type;  

private:

  std::unique_ptr<XCIntegratorImpl> pimpl_;

public:

  XCIntegraor( MPI_Comm, functional_type, const BasisSet&, const Molecule&,
    const MolGrid& );

};

}
