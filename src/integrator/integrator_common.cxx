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


//  for (int i = 0; i < submat_map.size(); i++) {
//      std::cout << i << "\t" << submat_map[i].first << "\t" << submat_map[i].second << std::endl;
//  }
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
  const int block_size = 256;
  std::vector< std::pair<int32_t, int32_t> > submat_map_expand;
  std::vector< int32_t > submat_block_idx;
  submat_block_idx.push_back(0);
  const int end_point = submat_map.back().second;

  int cut_index = 0;
  int cut_expand_index = 0;
  int small_index = 0;
  int delta;
  for (int block_start = 0; block_start < end_point; block_start += block_size) {
    const int block_end = block_start + block_size;
    
    int cut_start = submat_map[cut_index].first;
    int cut_end   = submat_map[cut_index].second;
    while (cut_index < submat_map.size() && cut_start < block_end) {
      if (cut_start < block_start && cut_end < block_start) {
	std::cout << "Something is wrong constructing the extended cut map " << cut_index << " " << cut_start << " "  << cut_end << " "  << block_start << std::endl;
      } else if (cut_start < block_start && cut_end > block_end) {
	submat_map_expand.emplace_back(block_start, block_end);

        delta = submat_map_expand.back().second - submat_map_expand.back().first;
        submat_map_expand.emplace_back(small_index, small_index + delta);
        small_index += delta;

	cut_expand_index++;
	break;
      } else if (cut_start < block_start) {
	submat_map_expand.emplace_back(block_start, cut_end);

        delta = submat_map_expand.back().second - submat_map_expand.back().first;
        submat_map_expand.emplace_back(small_index, small_index + delta);
        small_index += delta;

	cut_index++;
	cut_expand_index++;
      } else if (cut_end > block_end) {
	submat_map_expand.emplace_back(cut_start, block_end);

        delta = submat_map_expand.back().second - submat_map_expand.back().first;
        submat_map_expand.emplace_back(small_index, small_index + delta);
        small_index += delta;

	cut_expand_index++;
	break;
      } else {
        submat_map_expand.emplace_back(cut_start, cut_end);

        delta = submat_map_expand.back().second - submat_map_expand.back().first;
        submat_map_expand.emplace_back(small_index, small_index + delta);
        small_index += delta;

	cut_index++;
	cut_expand_index++;
      }

      cut_start = submat_map[cut_index].first;
      cut_end   = submat_map[cut_index].second;
    }
    submat_block_idx.push_back(cut_expand_index);
  }

  submat_map_expand[1].second = submat_block_idx[1];


  //for (int i = 0; i < submat_block_idx.size(); i++) {
  //  std::cout << submat_block_idx[i] << " ";
 // }
 // std::cout << std::endl;



/*
 // std::cout << "Task lengths: " << submat_map_check.size() << " " << submat_map_expand.size() << std::endl;
 // std::cout << "BlockIndex  : " << submat_map_check[1].second << " " << submat_map_expand[1].second << std::endl;
  for (int i = 0; i < submat_map_check.size(); i+=2) {
    if (submat_map_check[i].first != submat_map_expand[i].first)
      std::cout << "start:\t " << submat_map_check[i].first << " " << submat_map_expand[i].first << std::endl;
    if (submat_map_check[i].second != submat_map_expand[i].second)
      std::cout << "end  :\t " << submat_map_check[i].second << " " << submat_map_expand[i].second << std::endl;
    if (submat_map_check[i+1].first != submat_map_expand[i+1].first)
      std::cout << "small:\t " << submat_map_check[i+1].first << " " << submat_map_expand[i+1].first << std::endl;
  }

 // for (int i = submat_map_check.size(); i < submat_map_expand.size(); i+=2) {
 //     std::cout << "start:\t " << submat_map_expand[i].first << std::endl;
 // }

*/


  return submat_map_expand;
}



}
}
