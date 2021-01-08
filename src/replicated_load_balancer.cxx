#include "replicated_load_balancer.hpp"
#include <chrono>
#include <gauxc/util/div_ceil.hpp>

#ifdef GAUXC_ENABLE_CUDA
#include "load_balancer/cuda/cuda_collision_detection.hpp"
#endif

namespace GauXC {
namespace detail {

ReplicatedLoadBalancer::ReplicatedLoadBalancer( const ReplicatedLoadBalancer& ) = default;
ReplicatedLoadBalancer::ReplicatedLoadBalancer( ReplicatedLoadBalancer&& ) noexcept = default;

ReplicatedLoadBalancer::~ReplicatedLoadBalancer() noexcept = default;

std::unique_ptr<LoadBalancerImpl> ReplicatedLoadBalancer::clone() const {
  return std::make_unique<ReplicatedLoadBalancer>(*this);
}


auto raw_accumulate_tuple(
  const BasisSet<double>&      bs,
  size_t idx,
  const std::vector<int32_t>& counts,
  int32_t* raw_position_list
) {
  int32_t start = 0;
  if (idx != 0) start += counts[idx-1];
  int32_t end = counts[idx];

  std::vector<int32_t> shell_list(end - start);
  std::copy(raw_position_list + start, raw_position_list + end, shell_list.begin());

  size_t nbe = std::accumulate( shell_list.begin(), shell_list.end(), 0ul,
    [&](const auto& a, const auto& b) { return a + bs[b].size(); } );

  return std::tuple( std::move( shell_list ), nbe );
}


std::vector< XCTask > ReplicatedLoadBalancer::create_local_tasks_() const  {

  std::chrono::high_resolution_clock::time_point screen_start;
  std::chrono::high_resolution_clock::time_point screen_stop;
  double duration;

  const int32_t n_deriv = 1;

  int32_t world_rank;
  int32_t world_size;

#ifdef GAUXC_ENABLE_MPI
  MPI_Comm_rank( comm_, &world_rank );
  MPI_Comm_size( comm_, &world_size );
#else
  world_rank = 0;
  world_size = 1;
#endif

  std::vector< XCTask > local_work;
  std::vector<size_t> global_workload( world_size, 0 );   

  // Loop over Atoms
  const auto natoms = this->mol_->natoms();
  int32_t iCurrent  = 0;
  int32_t atBatchSz = 1;

  const size_t max_nbatches = mg_->max_nbatches();
  std::vector< std::pair<size_t, XCTask> > temp_tasks;
  temp_tasks.reserve( max_nbatches );

 // // Preallocate
 // std::vector<std::array<double,3>> batch_box_lo, batch_box_up;
 // if( screen_on_device ) {
 //   temp_tasks.resize(atBatchSz   * max_nbatches);
 //   batch_box_lo.resize(atBatchSz * max_nbatches);
 //   batch_box_up.resize(atBatchSz * max_nbatches);
 // }


  // Set up the collision problem
  std::vector<std::array<double,3>> low_points; 
  std::vector<std::array<double,3>> high_points;
  std::vector<std::array<double,3>> centers; centers.reserve((*this->basis_).size());
  std::vector<double> radii; radii.reserve((*this->basis_).size());

  screen_start = std::chrono::high_resolution_clock::now();
  for(auto& shell : (*this->basis_)) {
    centers.push_back(std::move(shell.O()));
    radii.push_back(shell.cutoff_radius());
  }
  for( const auto& atom : *this->mol_ ) {
    const std::array<double,3> center = { atom.x, atom.y, atom.z };

    auto& batcher = mg_->get_grid(atom.Z).batcher();
    batcher.quadrature().recenter( center );
    const size_t nbatches = batcher.nbatches();

    for( size_t ibatch = 0; ibatch < nbatches; ++ibatch ) {
      // Generate the batch (non-negligible cost)
      auto [lo, up, points, weights] = batcher.at(ibatch);

      low_points.push_back(std::move(lo));
      high_points.push_back(std::move(up));
    }
  }


  size_t LD_bit = util::div_ceil(centers.size(), 32);

#ifdef GAUXC_ENABLE_CUDA
  int32_t* raw_position_list;
  std::vector<int32_t> pos_list_idx(low_points.size());

  screen_start = std::chrono::high_resolution_clock::now();
  integrator::cuda::collision_detection(low_points.size(), centers.size(), LD_bit, low_points[0].data(), high_points[0].data(), centers[0].data(), radii.data(), pos_list_idx.data(), &raw_position_list);
  screen_stop = std::chrono::high_resolution_clock::now();
  duration = (double) std::chrono::duration_cast<std::chrono::milliseconds>( screen_stop - screen_start ).count();
  std::cout << "Generating cuda bitvector took: " << duration / 1000 << " seconds" << std::endl;
#endif


  // For batching of multiple atom screening
  size_t batch_idx_offset = 0;
  int idx = 0;

  screen_start = std::chrono::high_resolution_clock::now();
  for( const auto& atom : *this->mol_ ) {
    const std::array<double,3> center = { atom.x, atom.y, atom.z };

    auto& batcher = mg_->get_grid(atom.Z).batcher();
    batcher.quadrature().recenter( center );
    const size_t nbatches = batcher.nbatches();

    for( size_t ibatch = 0; ibatch < nbatches; ++ibatch ) {
    
      size_t batch_idx = ibatch + batch_idx_offset;

      // Generate the batch (non-negligible cost)
      auto [lo, up, points, weights] = batcher.at(ibatch);

      if( points.size() == 0 ) continue;

      auto [shell_list, nbe] = raw_accumulate_tuple( (*this->basis_), idx, pos_list_idx, raw_position_list);
      idx++;

      // Course grain screening
      if( not shell_list.size() ) continue; 

      // Copy task data
      XCTask task;
      task.iParent    = iCurrent;
      task.points     = std::move( points );
      task.weights    = std::move( weights );
      task.shell_list = std::move(shell_list);
      task.nbe        = nbe;
      task.dist_nearest = molmeta_->dist_nearest()[iCurrent];
      
      temp_tasks.push_back( 
        std::pair(batch_idx,std::move( task )) 
      );

    } // omp parallel for over batches





    // Assign Tasks to MPI ranks
    if( (iCurrent+1) % atBatchSz == 0 or iCurrent == ((int32_t)natoms-1) ) {

      // Sort based on task index for deterministic assignment
      std::sort( temp_tasks.begin(), temp_tasks.end(), 
        []( const auto& a, const auto& b ) {
          return a.first < b.first;
        } );

      // Assign batches to MPI ranks
      for( size_t ibatch = 0; ibatch < temp_tasks.size(); ++ibatch ) {

        XCTask task = std::move(temp_tasks.at(ibatch).second);
        //auto& points = task.points;
        //auto  nbe    = task.nbe;

        // Get rank with minimum work
        auto min_rank_it = 
          std::min_element( global_workload.begin(), global_workload.end() );
        int64_t min_rank = std::distance( global_workload.begin(), min_rank_it );

        global_workload[ min_rank ] += task.cost( n_deriv, natoms );

        if( world_rank == min_rank ) 
          local_work.push_back( std::move(task) );

      }

      temp_tasks.clear();

    }


    // Update counters and offsets
    iCurrent++;
    batch_idx_offset += nbatches;

  } // Loop over Atoms

  free(raw_position_list);

  screen_stop = std::chrono::high_resolution_clock::now();
  duration = (double) std::chrono::duration_cast<std::chrono::milliseconds>( screen_stop - screen_start ).count();
  std::cout << "Total loop: " << duration / 1000 << " seconds" << std::endl;

  screen_start = std::chrono::high_resolution_clock::now();
  // Lexicographic ordering of tasks
  auto task_order = []( const auto& a, const auto& b ) {

    // Sort by iParent first
    if( a.iParent < b.iParent )      return true;
    else if( a.iParent > b.iParent ) return false;

    // Equal iParent: lex sort on shell list
    else return a.shell_list < b.shell_list;

  };

  std::sort( local_work.begin(), local_work.end(),
    task_order ); 


  // Get unique tasks
  auto task_equiv = []( const auto& a, const auto& b ) {
    return a.equiv_with(b);
  };

  auto local_work_unique = local_work;
  auto last_unique = 
    std::unique( local_work_unique.begin(),
                 local_work_unique.end(),
                 task_equiv );
  local_work_unique.erase( last_unique, local_work_unique.end() );
  

  // Merge tasks
  for( auto&& t : local_work_unique ) {
    t.points.clear();
    t.weights.clear();
  }

  auto cur_lw_begin = local_work.begin();
  auto cur_uniq_it  = local_work_unique.begin();

  for( auto lw_it = local_work.begin(); lw_it != local_work.end(); ++lw_it ) 
  if( not task_equiv( *lw_it, *cur_uniq_it ) ) {

    if( cur_uniq_it == local_work_unique.end() )
      throw std::runtime_error("Messed up in unique");

    cur_uniq_it->merge_with( cur_lw_begin, lw_it );

    cur_lw_begin = lw_it;
    cur_uniq_it++;

  }

  // Merge the last set of batches
  for( ; cur_lw_begin != local_work.end(); ++cur_lw_begin )
    cur_uniq_it->merge_with( *cur_lw_begin );
  cur_uniq_it++;
  

  local_work = std::move(local_work_unique);
  
  screen_stop = std::chrono::high_resolution_clock::now();
  duration = (double) std::chrono::duration_cast<std::chrono::milliseconds>( screen_stop - screen_start ).count();
  std::cout << "Clean up: " << duration / 1000 << " seconds" << std::endl;
  return local_work;
}







}
}
