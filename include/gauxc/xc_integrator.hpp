#pragma once

#include <mpi.h>
#include <memory>

#include <gauxc/types.hpp>
#include <gauxc/basisset.hpp>
#include <gauxc/molecule.hpp>
#include <gauxc/molgrid.hpp>

namespace GauXC {

namespace detail {
  template <typename MatrixType>
  class XCIntegratorImpl;
}

template <typename MatrixType>
class XCIntegrator {

public:

  using matrix_type   = MatrixType;
  using value_type    = typename matrix_type::value_type;  
  using basisset_type = BasisSet< value_type >;

private:

  using pimpl_type    = detail::XCIntegratorImpl<MatrixType>;

  std::unique_ptr<pimpl_type> pimpl_;

public:

  XCIntegrator( pimpl_type&& pimpl );

  XCIntegrator( MPI_Comm, const functional_type&, const basisset_type&, 
    const Molecule&, const MolGrid& );

  XCIntegrator() noexcept;

};

}
