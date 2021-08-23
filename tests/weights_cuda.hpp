#pragma once
#include "weights_generate.hpp"
#include <fstream>
#include <string>
#include <gauxc/util/div_ceil.hpp>

#ifdef GAUXC_ENABLE_CUDA
#include <gauxc/util/cuda_util.hpp>
#include "device/cuda/cuda_aos_scheme1.hpp"
#include "device/cuda/cuda_aos_scheme1_weights.hpp"
      


                        
void test_cuda_weights( std::ifstream& in_file ) {

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

  constexpr auto weight_unroll = CudaAoSScheme1::weight_unroll;

  size_t LDatoms = util::div_ceil( natoms, weight_unroll ) * weight_unroll;

  std::vector< double >  coords( 3 * natoms );
  for( auto iat = 0 ; iat < natoms; ++iat ) {
    coords[ 3*iat + 0 ] = ref_data.mol.at(iat).x;
    coords[ 3*iat + 1 ] = ref_data.mol.at(iat).y;
    coords[ 3*iat + 2 ] = ref_data.mol.at(iat).z;
  }


  auto* points_d  = util::cuda_malloc<double>( 3*npts );
  auto* weights_d = util::cuda_malloc<double>( npts   );
  auto* iparent_d = util::cuda_malloc<int32_t>( npts  );
  auto* distnea_d = util::cuda_malloc<double>( npts   );
  auto* rab_d     = util::cuda_malloc<double>( natoms*LDatoms );
  auto* coords_d  = util::cuda_malloc<double>( 3*natoms );
  auto* dist_scr_d= util::cuda_malloc<double>( npts*LDatoms );

  util::cuda_copy( 3*npts, points_d,  points.data()->data() );
  util::cuda_copy( npts,   weights_d, weights.data() );
  util::cuda_copy( npts,   iparent_d, iparent.data() );
  util::cuda_copy( npts,   distnea_d, dist_nearest.data() );

  std::vector<double> rab_inv(natoms*natoms);
  for( auto i = 0; i < natoms*natoms; ++i )
    rab_inv[i] = 1./ref_data.meta->rab().data()[i];

  util::cuda_copy_2d( rab_d, LDatoms * sizeof(double),
                      rab_inv.data(), natoms * sizeof(double),
                      natoms * sizeof(double), natoms, "RAB H2D");

  util::cuda_copy( 3*natoms, coords_d, coords.data() );

  cudaStream_t stream = 0;
  cuda_aos_scheme1_weights_wrapper( npts, natoms, points_d, rab_d,
    LDatoms, coords_d, dist_scr_d, LDatoms, iparent_d, distnea_d,
    weights_d, stream );

  util::cuda_device_sync();
  util::cuda_copy( npts, weights.data(), weights_d );
  util::cuda_free( points_d, weights_d, iparent_d, distnea_d,
                   rab_d, coords_d, dist_scr_d );

  for( auto i = 0ul; i < npts; ++i )
    CHECK( weights.at(i) == Approx( weights_ref.at(i) ) );

}
#endif
