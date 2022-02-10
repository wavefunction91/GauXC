#pragma once

#include <gauxc/molgrid.hpp>

namespace GauXC {
namespace detail {

class MolGridImpl {

  atomic_grid_map molgrid_;

public:

  MolGridImpl( const atomic_grid_map& );

  MolGridImpl( const MolGridImpl& );
  MolGridImpl( MolGridImpl&& ) noexcept;

  ~MolGridImpl() noexcept;

  size_t natoms_uniq() const;
  const Grid& get_grid( AtomicNumber ) const;
        Grid& get_grid( AtomicNumber )      ;

  size_t max_nbatches() const;

};

}
}
