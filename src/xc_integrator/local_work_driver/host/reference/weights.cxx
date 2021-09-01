#include "host/reference/weights.hpp"
#include "common/integrator_constants.hpp"

namespace GauXC {

// Reference Becke weights impl
void reference_becke_weights_host(
  const Molecule&        mol,
  const MolMeta&         meta,
  task_iterator          task_begin,
  task_iterator          task_end
) {

  // Becke partition functions
  auto hBecke = [](double x) {return 1.5 * x - 0.5 * x * x * x;}; // Eq. 19
  auto gBecke = [&](double x) {return hBecke(hBecke(hBecke(x)));}; // Eq. 20 f_3

  const size_t ntasks = std::distance(task_begin,task_end);
  const size_t natoms = mol.natoms();

  const auto&  RAB    = meta.rab();

  #pragma omp parallel 
  {

  std::vector<double> partitionScratch( natoms );
  std::vector<double> atomDist( natoms );

  #pragma omp for
  for( size_t iT = 0; iT < ntasks;                  ++iT )
  for( size_t i  = 0; i  < (task_begin+iT)->points.size(); ++i  ) {

    auto&       task   = *(task_begin+iT);
    auto&       weight = task.weights[i];
    const auto& point  = task.points[i];

    // Compute distances of each center to point
    for(size_t iA = 0; iA < natoms; iA++) {

      const double da_x = point[0] - mol[iA].x;
      const double da_y = point[1] - mol[iA].y;
      const double da_z = point[2] - mol[iA].z;

      atomDist[iA] = std::sqrt(da_x*da_x + da_y*da_y + da_z*da_z);

    }

    // Evaluate unnormalized partition functions 
    std::fill(partitionScratch.begin(),partitionScratch.end(),1.);
    for( size_t iA = 0; iA < natoms; iA++ ) 
    for( size_t jA = 0; jA < iA;     jA++ ){
      const double mu = (atomDist[iA] - atomDist[jA]) / RAB[jA + iA*natoms];
      const double g = gBecke(mu);

      partitionScratch[iA] *= 0.5 * (1. - g);
      partitionScratch[jA] *= 0.5 * (1. + g);
    }

    // Normalization
    double sum = 0.;
    for( size_t iA = 0; iA < natoms; iA++ )  sum += partitionScratch[iA];

    // Update Weights
    weight *= partitionScratch[task.iParent] / sum;

  } // Collapsed loop over tasks and points

  } // OMP context


}

void reference_ssf_weights_host(
  const Molecule&        mol,
  const MolMeta&         meta,
  task_iterator          task_begin,
  task_iterator          task_end
) {

  auto gFrisch = [&](double x) {
    const double s_x  = x / integrator::magic_ssf_factor<>;
    const double s_x2 = s_x  * s_x;
    const double s_x3 = s_x  * s_x2;
    const double s_x5 = s_x3 * s_x2;
    const double s_x7 = s_x5 * s_x2;

    return (35.*(s_x - s_x3) + 21.*s_x5 - 5.*s_x7) / 16.;
  };

  const size_t ntasks = std::distance(task_begin,task_end);
  const size_t natoms = mol.natoms();

  const auto&  RAB    = meta.rab();

  #pragma omp parallel 
  {

  std::vector<double> partitionScratch( natoms );
  std::vector<double> atomDist( natoms );

  #pragma omp for
  for( size_t iT = 0; iT < ntasks;                  ++iT )
  for( size_t i  = 0; i  < (task_begin+iT)->points.size(); ++i  ) {

    auto&       task   = *(task_begin+iT);
    auto&       weight = task.weights[i];
    const auto& point  = task.points[i];

    const auto dist_cutoff = 0.5 * (1-integrator::magic_ssf_factor<>) * task.dist_nearest;

    // Compute dist to parent atom
    {
      const double da_x = point[0] - mol[task.iParent].x;
      const double da_y = point[1] - mol[task.iParent].y;
      const double da_z = point[2] - mol[task.iParent].z;

      atomDist[task.iParent] = std::sqrt(da_x*da_x + da_y*da_y + da_z*da_z);
    }

    if( atomDist[task.iParent] < dist_cutoff ) continue; // Partition weight = 1

    // Compute distances of each center to point
    for(size_t iA = 0; iA < natoms; iA++) {

      if( iA == (size_t)task.iParent ) continue;

      const double da_x = point[0] - mol[iA].x;
      const double da_y = point[1] - mol[iA].y;
      const double da_z = point[2] - mol[iA].z;

      atomDist[iA] = std::sqrt(da_x*da_x + da_y*da_y + da_z*da_z);

    }

    // Evaluate unnormalized partition functions 
    std::fill(partitionScratch.begin(),partitionScratch.end(),1.);
    for( size_t iA = 0; iA < natoms; iA++ ) 
    for( size_t jA = 0; jA < iA;     jA++ )
    if( partitionScratch[iA] > integrator::ssf_weight_tol or 
        partitionScratch[jA] > integrator::ssf_weight_tol ) {

      const double mu = (atomDist[iA] - atomDist[jA]) / RAB[jA + iA*natoms];

      if( mu <= -integrator::magic_ssf_factor<> ) {

        partitionScratch[jA] = 0.;

      } else if (mu >= integrator::magic_ssf_factor<>) {

        partitionScratch[iA] = 0.;

      } else {

        double g = 0.5 * ( 1. - gFrisch(mu) );
        partitionScratch[iA] *= g;
        partitionScratch[jA] *= 1. - g;

      }

    }

    // Normalization
    double sum = 0.;
    for( size_t iA = 0; iA < natoms; iA++ )  sum += partitionScratch[iA];

    // Update Weights
    weight *= partitionScratch[task.iParent] / sum;

  } // Collapsed loop over tasks and points

  } // OMP context


}

void reference_lko_weights_host(
  const Molecule&        mol,
  const MolMeta&         meta,
  task_iterator          task_begin,
  task_iterator          task_end
) {

#if 0
  // Sort on atom index
  std::stable_sort( task_begin, task_end, [](const auto& a, const auto&b ) { return a.iParent < b.iParent; } );

  auto RAB = meta.rab();
  const auto natoms = mol.size();

  auto it = task_begin;
  for( auto iAtom = 0; iAtom < natoms; ++iAtom ) {

    auto atom_end = std::find_if( it, task_end, [&](const auto& a){ return a.iParent == (iAtom+1); } );
    
    size_t npts_atom = std::accumulate( it, atom_end, 0ul, [&](const auto& a, const auto&b) {
      return a + b.points.size();
    });

    std::vector<size_t> atom_idx(natoms);
    std::iota( atom_idx.begin(), atom_idx.end(), 0 );
    std::sort( atom_idx.begin(), atom_idx.end(), [&](auto i, auto j) {
      return RAB[i + iAtom*natoms] < RAB[j + iAtom*natoms];
    });
    if( atom_idx[0] != iAtom ) throw std::runtime_error("WTF");

    std::cout << "Parent Atom " << iAtom << std::endl;
    for( auto i = 0; i < natoms; ++i ) {
      std::cout << "  " << atom_idx[i] << ", " << RAB[atom_idx[i] + iAtom*natoms] << std::endl;
    }
    std::cout << std::endl;

    constexpr double R_cutoff = 5;
    auto h_func = [](auto x){ return 1.5*x - 0.5 * x*x*x; };
    auto g_func = [&](auto x){ return h_func(h_func(h_func(x))); };

    // Loop over tasks / weights / points
    for( auto task_it = it; task_it != atom_end; task_it++ ) {
      auto& points  = task_it->points;
      auto& weights = task_it->weights;
      const auto npts = points.size();
    for( auto ipt = 0; ipt < npts; ++ipt ) {

      

      // Compute distance from point to parent atom
      const double r_parent_x = mol[iAtom].x - points[ipt][0];
      const double r_parent_y = mol[iAtom].y - points[ipt][1];
      const double r_parent_z = mol[iAtom].z - points[ipt][2];

      const auto r_parent = std::sqrt(
        r_parent_x*r_parent_x + r_parent_y*r_parent_y + r_parent_z*r_parent_z
      );

      double r_nearest = r_parent;
      std::vector<double> r_atoms(natoms, std::numeric_limits<double>::max()); r_atoms[iAtom] = r_parent;
      size_t natoms_keep = 1;
      for( auto _idx = 1; _idx < natoms; _idx++ ) {
        auto idx = atom_idx[_idx];
        auto Rpi = RAB[ idx + iAtom*natoms ];
        if( Rpi > r_nearest + r_parent + 2*R_cutoff ) break;

        // Compute distance from point to this atom
        const double r_atoms_x = mol[idx].x - points[ipt][0];
        const double r_atoms_y = mol[idx].y - points[ipt][1];
        const double r_atoms_z = mol[idx].z - points[ipt][2];

        r_atoms[idx] = std::sqrt(
          r_atoms_x*r_atoms_x + r_atoms_y*r_atoms_y + r_atoms_z*r_atoms_z
        );

        r_nearest = std::min( r_nearest, r_atoms[idx] );
        natoms_keep++;
      }

      std::cout << "  R PARENT  = " << r_parent << std::endl;
      std::cout << "  R NEAREST = " << r_nearest << std::endl;

      //if( r_nearest > r_parent + R_cutoff ) continue; // XXX This could be wrong, it's f-ed up in the paper...
      //if( r_parent > r_nearest + R_cutoff ) continue; // XXX This could be wrong, it's f-ed up in the paper...
      if( r_nearest > r_parent + R_cutoff ) {
        weights[ipt] = 0.;
        continue;
      }


      std::vector<size_t> point_idx(natoms);
      std::iota(point_idx.begin(),point_idx.end(),0);
      std::sort( point_idx.begin(), point_idx.end(), [&](auto i, auto j) {
        return r_atoms[i] < r_atoms[j];
      } );

      std::vector<double> P(natoms,-1.);
      for( auto i_idx = 0; i_idx < natoms; ++i_idx ) {
        auto pt_i = point_idx[i_idx];
        auto at_i = atom_idx[pt_i];
        if( r_atoms[pt_i] > r_nearest + R_cutoff ) break;
        for( auto j_idx = i_idx+1; j_idx < natoms; ++j_idx ) {
          auto pt_j = point_idx[j_idx];
          auto at_j = atom_idx[pt_j];
          if( r_atoms[pt_j] > r_atoms[pt_i] + R_cutoff ) break;

          const auto mu = (r_atoms[pt_i] - r_atoms[pt_j]) / RAB[ at_i + at_j*natoms];
          auto g_mu = g_func(mu);
          auto s_ij = 0.5 * ( 1. - g_mu );
          auto s_ji = 1. - s_ij;
          if(P[pt_i] < 0.) P[pt_i] =  s_ij;
          else             P[pt_i] *= s_ij;
          if(P[pt_j] < 0.) P[pt_j] =  s_ji;
          else             P[pt_j] *= s_ji;
        }
      }

      for( auto& p : P ) if( p < 0. ) p = 0.;

      double sum = std::accumulate( P.begin(), P.end(), 0. );
      weights[ipt] *= P[iAtom] / sum;


    }
    }

    it = atom_end;

  }
#else

  // Sort on atom index
  std::stable_sort( task_begin, task_end, 
    [](const auto& a, const auto&b ) { return a.iParent < b.iParent; } );

  // Becke partition functions
  auto hBecke = [](double x) {return 1.5 * x - 0.5 * x * x * x;}; // Eq. 19
  auto gBecke = [&](double x) {return hBecke(hBecke(hBecke(x)));}; // Eq. 20 f_3
  constexpr double R_cutoff = 5;

  const size_t natoms = mol.natoms();

  const auto&  RAB    = meta.rab();

  //#pragma omp parallel 
  {

  std::vector<double> partitionScratch( natoms );
  std::vector<double> atomDist( natoms );
  std::vector<size_t> inter_atom_dist_idx( natoms );
  std::vector<size_t> point_dist_idx( natoms );

  //#pragma omp for schedule(dynamic)
  for( auto iAtom = 0; iAtom < natoms; ++iAtom ) {

    auto atom_begin = std::find_if( task_begin, task_end,
      [&](const auto& t){ return t.iParent == iAtom; } );
    auto atom_end = std::find_if( task_begin, task_end,
      [&](const auto& t){ return t.iParent == (iAtom+1); } );

    auto* RAB_parent = RAB.data() + iAtom*natoms;

    std::iota( inter_atom_dist_idx.begin(), inter_atom_dist_idx.end(), 0 );
    std::sort( inter_atom_dist_idx.begin(), inter_atom_dist_idx.end(),
      [&](auto i, auto j){ return RAB_parent[i] < RAB_parent[j]; } );

  for( auto task_it = atom_begin; task_it != atom_end; ++task_it ) {

    auto& points  = task_it->points;
    auto& weights = task_it->weights;
    const auto npts = points.size();

  for( auto ipt = 0; ipt < npts; ++ipt ) {

    auto& weight = weights[ipt];
    const auto point = points[ipt];

    std::fill( atomDist.begin(), atomDist.end(), std::numeric_limits<double>::infinity() );
    // Parent distance
    {
      const double da_x = point[0] - mol[iAtom].x;
      const double da_y = point[1] - mol[iAtom].y;
      const double da_z = point[2] - mol[iAtom].z;

      atomDist[iAtom] = std::sqrt(da_x*da_x + da_y*da_y + da_z*da_z);
    }

    double r_parent  = atomDist[iAtom];
    double r_nearest = r_parent;
    // Compute distances of each center to point
    for(size_t iA = 1; iA < natoms; iA++) {
      auto idx = inter_atom_dist_idx[iA];
      if( RAB_parent[idx] > (r_parent + r_nearest + 2*R_cutoff) ) break;

      const double da_x = point[0] - mol[idx].x;
      const double da_y = point[1] - mol[idx].y;
      const double da_z = point[2] - mol[idx].z;

      atomDist[idx] = std::sqrt(da_x*da_x + da_y*da_y + da_z*da_z);
      r_nearest = std::min( r_nearest, atomDist[idx] );
    }

     // Partition weight is 0
    if( r_parent > r_nearest + R_cutoff ) {
      weight = 0.;
      continue;
    }

    std::iota( point_dist_idx.begin(), point_dist_idx.end(), 0 );
    std::sort( point_dist_idx.begin(), point_dist_idx.end(),
      [&](auto i, auto j){ return atomDist[i] < atomDist[j]; } );



    // Evaluate unnormalized partition functions 
    std::fill(partitionScratch.begin(),partitionScratch.end(),0.);
    for( int i = 0; i < natoms; ++i ) {
      auto idx_i = point_dist_idx[i];
      auto r_i = atomDist[idx_i];
      if( r_i > (r_nearest + R_cutoff) ) { break; }
      partitionScratch[idx_i] = 1.;
    for( int j = 0; j < i; ++j ) {
      auto idx_j = point_dist_idx[j];
      auto r_j = atomDist[idx_j];
      if( r_j > (r_i + R_cutoff) ) { break; }

      const double mu = 
        (r_i - r_j) / std::min(RAB[idx_j + idx_i*natoms], R_cutoff);
      if( mu <= -1 ) {
        partitionScratch[idx_j] = 0.;
        continue;
      }

      if( mu >= 1 ) {
        partitionScratch[idx_i] = 0.;
        continue;
      }

      const double g = gBecke(mu);
      partitionScratch[idx_i] *= 0.5 * (1. - g);
      partitionScratch[idx_j] *= 0.5 * (1. + g);
      
    }
    }

    // Normalization
    double sum = 0.;
    for( size_t iA = 0; iA < natoms; iA++ )  sum += partitionScratch[iA];

    // Update Weights
    weight *= partitionScratch[iAtom] / sum;

  } // Loop over points 
  } // Loop over tasks
  } // Loop over atoms

  } // OMP context

#endif
}

}
