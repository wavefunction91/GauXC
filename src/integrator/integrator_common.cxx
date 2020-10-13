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


  /*
   * This code block does post-processing for the submatrix optimizations
   *
   * It first adds the index within the small matrix as another pair in the vector.
   * This allows the kernel to process multiple cuts concurrently within the same
   * task. Additionally, it adds artificial breaks in the cut at the given interval
   * This is to reduce the amount of bookkeeping that the kernel is required to do.
   *
   * While the small matrix start indices are stored in the additional pair, the second 
   * value is blank as the delta can be reused from the big matrix start and stop points.
   *
   * The L2 cache optimization forces breaks this into muliple kernel calls, and I am 
   * currently using the forth value in a cut to pass information about the block to the kernel.
   * This is veyr temporary and primarily so I did not have to add an additional input variable.
   *
   */
  const int block_size = 512;
  std::vector< std::pair<int32_t, int32_t> > submat_map_expand;
  int small_index = 0;
  bool passedBlock = false;
  int delta;
  for (int i = 0; i < submat_map.size(); i++) {
    int start = submat_map[i].first;
    int end   = submat_map[i].second;
    while ((start / block_size) < (end / block_size)) {
      int next_end = ((start / block_size) + 1) * block_size;

      submat_map_expand.emplace_back(start, next_end);

      delta = next_end - start;
      submat_map_expand.emplace_back(small_index, small_index + delta);
      small_index += delta;

      if (start >= block_size && !passedBlock) {
        submat_map_expand[1].second = i;
	passedBlock = true;
      }

      start = next_end;
    }

    if (start != end) {
       submat_map_expand.emplace_back(start, end);

      delta = end - start;
      submat_map_expand.emplace_back(small_index, small_index + delta);
      small_index += delta;


      if (start >= block_size && !passedBlock) {
        submat_map_expand[1].second = i;
	passedBlock = true;
      }
    }

  }


  return submat_map_expand;
}



}
}
