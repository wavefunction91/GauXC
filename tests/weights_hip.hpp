#pragma once
#include "weights_generate.hpp"
#include <fstream>
#include <string>
#include <gauxc/util/div_ceil.hpp>

#ifdef GAUXC_ENABLE_HIP
#include <gauxc/util/hip_util.hpp>
#include "device/hip/kernels/grid_to_center.hpp"
#include "device/hip/kernels/hip_ssf_1d.hpp"
      


                        
void test_hip_weights( std::ifstream& in_file ) {

  ref_weights_data ref_data;
  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }

  std::vector< std::array<double,3> > points;
  std::vector< double >               weights, weights_ref;
  std::vector< double >               dist_nearest;
  std::vector< int32_t >              iparent;

  for( auto& task : ref_data.tasks_unm ) {
    points.insert( points.end(),
                   task.points.begin(),
                   task.points.end() );
    weights.insert( weights.end(),
                    task.weights.begin(),
                    task.weights.end() );

    size_t npts = task.points.size();
    dist_nearest.insert( dist_nearest.end(), npts,
                         task.dist_nearest );
    iparent.insert( iparent.end(), npts, task.iParent );
  }

  for( auto& task : ref_data.tasks_mod ) {
    weights_ref.insert( weights_ref.end(),
                        task.weights.begin(),
                        task.weights.end() );
  }

  size_t npts   = points.size();
  size_t natoms = ref_data.mol.natoms();

  size_t LDatoms = natoms;

  std::vector< double >  coords( 3 * natoms );
  for( auto iat = 0 ; iat < natoms; ++iat ) {
    coords[ 3*iat + 0 ] = ref_data.mol.at(iat).x;
    coords[ 3*iat + 1 ] = ref_data.mol.at(iat).y;
    coords[ 3*iat + 2 ] = ref_data.mol.at(iat).z;
  }


  auto* points_d  = util::hip_malloc<double>( 3*npts );
  auto* weights_d = util::hip_malloc<double>( npts   );
  auto* iparent_d = util::hip_malloc<int32_t>( npts  );
  auto* distnea_d = util::hip_malloc<double>( npts   );
  auto* rab_d     = util::hip_malloc<double>( natoms*LDatoms );
  auto* coords_d  = util::hip_malloc<double>( 3*natoms );
  auto* dist_scr_d= util::hip_malloc<double>( npts*LDatoms );

  util::hip_copy( 3*npts, points_d,  points.data()->data() );
  util::hip_copy( npts,   weights_d, weights.data() );
  util::hip_copy( npts,   iparent_d, iparent.data() );
  util::hip_copy( npts,   distnea_d, dist_nearest.data() );

  std::vector<double> rab_inv(natoms*natoms);
  for( auto i = 0; i < natoms*natoms; ++i )
    rab_inv[i] = 1./ref_data.meta->rab().data()[i];

  util::hip_copy_2d( rab_d, LDatoms * sizeof(double),
                      rab_inv.data(), natoms * sizeof(double),
                      natoms * sizeof(double), natoms, "RAB H2D");

  util::hip_copy( 3*natoms, coords_d, coords.data() );

  hipStream_t stream = 0;
#if 0
  hip_aos_scheme1_weights_wrapper( npts, natoms, points_d, rab_d,
    LDatoms, coords_d, dist_scr_d, LDatoms, iparent_d, distnea_d,
    weights_d, stream );
#else
  compute_grid_to_center_dist( npts, natoms, coords_d, points_d, dist_scr_d,
     LDatoms, stream );
  partition_weights_ssf_1d( npts, natoms, rab_d, LDatoms, coords_d, dist_scr_d,
    LDatoms, iparent_d, distnea_d, weights_d, stream );
#endif

  util::hip_device_sync();
  util::hip_copy( npts, weights.data(), weights_d );

  std::vector<double> dist_scr_host(natoms * npts);
  util::hip_copy( npts*natoms, dist_scr_host.data(), dist_scr_d ); 
  for( int i = 0; i < npts; ++i ) {
  for( int j = 0; j < natoms; ++j ) {

    const double x = points[i][0];
    const double y = points[i][1];
    const double z = points[i][2];

    const double Rx = coords[3*j + 0];
    const double Ry = coords[3*j + 1];
    const double Rz = coords[3*j + 2];

    const double rx = Rx - x;
    const double ry = Ry - y;
    const double rz = Rz - z;

    const double r = std::sqrt( rx*rx + ry*ry + rz*rz );
    CHECK( r == Approx(dist_scr_host[ j + i*LDatoms ]) );
  }
  }

  util::hip_free( points_d, weights_d, iparent_d, distnea_d,
                   rab_d, coords_d, dist_scr_d );

  for( auto i = 0ul; i < npts; ++i )
    CHECK( weights.at(i) == Approx( weights_ref.at(i) ) );

}
#endif
