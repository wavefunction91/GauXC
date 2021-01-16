#include "replicated_load_balancer.hpp"

namespace GauXC {
namespace detail {

HostReplicatedLoadBalancer::HostReplicatedLoadBalancer( const HostReplicatedLoadBalancer& ) = default;
HostReplicatedLoadBalancer::HostReplicatedLoadBalancer( HostReplicatedLoadBalancer&& ) noexcept = default;

HostReplicatedLoadBalancer::~HostReplicatedLoadBalancer() noexcept = default;

std::unique_ptr<LoadBalancerImpl> HostReplicatedLoadBalancer::clone() const {
  return std::make_unique<HostReplicatedLoadBalancer>(*this);
}

bool cube_sphere_intersect( 
  const std::array<double,3>& lo, 
  const std::array<double,3>& up,
  const std::array<double,3>& center,
  const double                rad
) {

  double dist = rad * rad;

  if( center[0] < lo[0] ) {
    const double r_lo = center[0] - lo[0];
    const double dist_lo = r_lo * r_lo;
    dist -= dist_lo;
  } else if( center[0] > up[0] ) {
    const double r_up = center[0] - up[0];
    const double dist_up = r_up * r_up;
    dist -= dist_up;
  }

  if( dist < 0. ) return false;

  if( center[1] < lo[1] ) {
    const double r_lo = center[1] - lo[1];
    const double dist_lo = r_lo * r_lo;
    dist -= dist_lo;
  } else if( center[1] > up[1] ) {
    const double r_up = center[1] - up[1];
    const double dist_up = r_up * r_up;
    dist -= dist_up;
  }

  if( dist < 0. ) return false;


  if( center[2] < lo[2] ) {
    const double r_lo = center[2] - lo[2];
    const double dist_lo = r_lo * r_lo;
    dist -= dist_lo;
  } else if( center[2] > up[2] ) {
    const double r_up = center[2] - up[2];
    const double dist_up = r_up * r_up;
    dist -= dist_up;
  }

  return dist > 0.;

}


auto micro_batch_screen(
  const BasisSet<double>&      bs,
  const std::array<double,3>&  box_lo,
  const std::array<double,3>&  box_up
) {


  std::vector<int32_t> shell_list; shell_list.reserve(bs.nshells());
  for(auto iSh = 0ul; iSh < bs.size(); ++iSh) {

    const auto& center = bs[iSh].O();
    const auto  crad   = bs[iSh].cutoff_radius();
    const bool intersect = 
      cube_sphere_intersect( box_lo, box_up, center, crad );
    

    //std::cout << "  MBS: " << iSh << ", " << 
    //          center[0] << ", " << center[1] << ", " << center[2] << ", " <<
    //          box_up[0] << ", " << box_up[1] << ", " << box_up[2] << ", " <<
    //          box_lo[0] << ", " << box_lo[1] << ", " << box_lo[2] << ", " <<
    //          crad << std::boolalpha << ", " << intersect << std::endl;
              

    // Add shell to list if need be
    if( intersect )
      shell_list.emplace_back( iSh );

  }

  size_t nbe = std::accumulate( shell_list.begin(), shell_list.end(), 0ul,
    [&](const auto& a, const auto& b) { return a + bs[b].size(); } );

  return std::tuple( std::move( shell_list ), nbe );

}
  






std::vector< XCTask > HostReplicatedLoadBalancer::create_local_tasks_() const  {

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

  // For batching of multiple atom screening
  size_t batch_idx_offset = 0;

  for( const auto& atom : *this->mol_ ) {


    const std::array<double,3> center = { atom.x, atom.y, atom.z };

    auto& batcher = mg_->get_grid(atom.Z).batcher();
    batcher.quadrature().recenter( center );




    const size_t nbatches = batcher.nbatches();

    #pragma omp parallel for
    for( size_t ibatch = 0; ibatch < nbatches; ++ibatch ) {
    
      size_t batch_idx = ibatch + batch_idx_offset;

      // Generate the batch (non-negligible cost)
      auto [lo, up, points, weights] = batcher.at(ibatch);

      if( points.size() == 0 ) continue;

      // Microbatch Screening
      auto [shell_list, nbe] = micro_batch_screen( (*this->basis_), lo, up );

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

      #pragma omp critical
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

  return local_work;
}







}
}
