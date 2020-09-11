#include "integrator_common.hpp"


namespace GauXC      {
namespace integrator {

std::vector< std::pair<int32_t, int32_t> >
  gen_compressed_submat_map( const BasisSet<double>&       basis,
                             const std::vector< int32_t >& shell_mask ) {


  std::vector< std::pair<int32_t, int32_t> > submat_map;

  // Init as if there is no screening
  submat_map.emplace_back(
    basis.shell_to_ao_range( shell_mask.front() ).first,
    basis.shell_to_ao_range( shell_mask.back()  ).second
  );


  for( auto sh_it =  shell_mask.begin(); sh_it != shell_mask.end()-1; ++sh_it ) {

    if( *(sh_it+1) - *(sh_it) != 1 ) {

      submat_map.back().second = basis.shell_to_ao_range(*sh_it).second;
        
      submat_map.emplace_back(
        basis.shell_to_ao_range( *(sh_it+1) ).first,
        basis.shell_to_ao_range( shell_mask.back()  ).second
      );

    }



  }


  if( shell_mask.size() == 1 )
    submat_map.back().second = 
      basis.shell_to_ao_range(shell_mask[0]).second;

  return submat_map;

}



}
}
