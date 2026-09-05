/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy).
 *
 * (c) 2024-2025, Microsoft Corporation
 *
 * All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include "weights_generate.hpp"
#include "hdf5_test_serialization.hpp"
#include "hdf5_test_serialization_impl.hpp"
#include <string>
#include <gauxc/util/div_ceil.hpp>

#ifdef GAUXC_HAS_SYCL
#include "device_specific/sycl_util.hpp"
#include "device/sycl/sycl_aos_scheme1.hpp"
#include "device/sycl/sycl_aos_scheme1_weights.hpp"




void test_sycl_weights( const std::string& filename ) {

  ref_weights_data ref_data;
  read_weights_data(ref_data, filename);

  //std::vector< std::array<double,3> > points;
  std::vector< double > points_x, points_y, points_z;
  std::vector< double >               weights, weights_ref;
  std::vector< double >               dist_nearest;
  std::vector< int32_t >              iparent;

  for( auto& task : ref_data.tasks_unm ) {
    for( auto pt : task.points ) {
      points_x.emplace_back(pt[0]);
      points_y.emplace_back(pt[1]);
      points_z.emplace_back(pt[2]);
    }
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

  size_t npts   = points_x.size();
  size_t natoms = ref_data.mol.natoms();

  constexpr auto weight_unroll = alg_constants::SyclAoSScheme1::weight_unroll;

  size_t LDatoms = util::div_ceil( natoms, weight_unroll ) * weight_unroll;

  std::vector< double >  coords( 3 * natoms );
  for( auto iat = 0 ; iat < natoms; ++iat ) {
    coords[ 3*iat + 0 ] = ref_data.mol.at(iat).x;
    coords[ 3*iat + 1 ] = ref_data.mol.at(iat).y;
    coords[ 3*iat + 2 ] = ref_data.mol.at(iat).z;
  }

  // SYCL USM allocation is context-bound, so every alloc/copy/free below
  // needs the queue that owns the allocations.
  util::sycl_queue stream;
  ::sycl::queue& q = stream.queue;

  //auto* points_d  = util::sycl_malloc<double>( 3*npts, q );
  auto* points_x_d  = util::sycl_malloc<double>( npts, q );
  auto* points_y_d  = util::sycl_malloc<double>( npts, q );
  auto* points_z_d  = util::sycl_malloc<double>( npts, q );
  auto* weights_d = util::sycl_malloc<double>( npts   , q );
  auto* iparent_d = util::sycl_malloc<int32_t>( npts  , q );
  auto* distnea_d = util::sycl_malloc<double>( npts   , q );
  auto* rab_d     = util::sycl_malloc<double>( natoms*LDatoms, q );
  auto* coords_d  = util::sycl_malloc<double>( 3*natoms, q );
  auto* dist_scr_d= util::sycl_malloc<double>( npts*LDatoms, q );

  //util::sycl_copy( 3*npts, points_d,  points.data()->data(), q );
  util::sycl_copy( npts, points_x_d,  points_x.data(), q );
  util::sycl_copy( npts, points_y_d,  points_y.data(), q );
  util::sycl_copy( npts, points_z_d,  points_z.data(), q );
  util::sycl_copy( npts,   weights_d, weights.data(), q );
  util::sycl_copy( npts,   iparent_d, iparent.data(), q );
  util::sycl_copy( npts,   distnea_d, dist_nearest.data(), q );

  std::vector<double> rab_inv(natoms*natoms);
  for( auto i = 0; i < natoms*natoms; ++i )
    rab_inv[i] = 1./ref_data.meta->rab().data()[i];

  util::sycl_copy_2d( rab_d, LDatoms * sizeof(double),
                      rab_inv.data(), natoms * sizeof(double),
                      natoms * sizeof(double), natoms, q, "RAB H2D");

  util::sycl_copy( 3*natoms, coords_d, coords.data(), q );

  sycl_aos_scheme1_weights_wrapper( npts, natoms, points_x_d, points_y_d, points_z_d, rab_d,
    LDatoms, coords_d, dist_scr_d, LDatoms, iparent_d, distnea_d,
    weights_d, q );

  util::sycl_device_sync(q);
  util::sycl_copy( npts, weights.data(), weights_d, q );
  util::sycl_free( q, points_x_d, points_y_d, points_z_d, weights_d, iparent_d, distnea_d,
                   rab_d, coords_d, dist_scr_d );

  for( auto i = 0ul; i < npts; ++i )
    CHECK( weights.at(i) == Approx( weights_ref.at(i) ) );

}
#endif
