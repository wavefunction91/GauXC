#pragma once
#include <gauxc/gauxc_config.hpp>
#include <gauxc/xc_integrator/xc_integrator_impl.hpp>
#include <gauxc/xc_integrator/xc_cuda_data.hpp>
#include <gauxc/xc_integrator/xc_cuda_util.hpp>

#ifdef GAUXC_ENABLE_CUDA
namespace GauXC  {
namespace detail {

using namespace GauXC::integrator::cuda;


template <typename MatrixType>
class ShellBatchedXCCudaIntegrator : public XCIntegratorImpl<MatrixType> {

  using base_type     = XCIntegratorImpl<MatrixType>;
  using matrix_type   = typename base_type::matrix_type;
  using value_type    = typename base_type::value_type;
  using basisset_type = typename base_type::basisset_type;
  using exc_vxc_type  = typename base_type::exc_vxc_type;
    
  std::shared_ptr< XCCudaData< value_type > > cuda_data_;

  exc_vxc_type eval_exc_vxc_( const MatrixType& ) override; 

public:

  template <typename... Args>
  ShellBatchedXCCudaIntegrator( Args&&... args ) :
    base_type( std::forward<Args>(args)... ) { }

  ShellBatchedXCCudaIntegrator( const ShellBatchedXCCudaIntegrator& ) = default;
  ShellBatchedXCCudaIntegrator( ShellBatchedXCCudaIntegrator&& ) noexcept = default;

  ~ShellBatchedXCCudaIntegrator() noexcept = default;

};




template <typename MatrixType>
typename ShellBatchedXCCudaIntegrator<MatrixType>::exc_vxc_type 
  ShellBatchedXCCudaIntegrator<MatrixType>::eval_exc_vxc_( const MatrixType& P ) {


#ifdef GAUXC_ENABLE_MAGMA
  // Initialize MAGMA
  {
    auto ierr = magma_init();
    GAUXC_MAGMA_ERROR( "MAGMA Init Failed", ierr );
  }
#endif

#ifdef GAUXC_ENABLE_MPI
  int32_t device_count, cur_device;
  cudaGetDeviceCount( &device_count );
  cudaGetDevice( &cur_device );
 
  int32_t world_rank, world_size;
  MPI_Comm_rank( this->comm_, &world_rank );
  MPI_Comm_size( this->comm_, &world_size );

/* XXX: Does not work on Summit
  MPI_Comm node_comm;
  MPI_Comm_split_type(this->comm_, MPI_COMM_TYPE_SHARED, 0,
                      MPI_INFO_NULL, &node_comm);

  int32_t node_rank, node_size;
  MPI_Comm_rank( node_comm, &node_rank );
  MPI_Comm_size( node_comm, &node_size );

  if( node_size > device_count )
    throw std::runtime_error("GauXC + CUDA Assumes MPI <-> GPU is 1-to-1");

  cudaSetDevice( node_rank );
*/
#endif


  size_t nbf     = this->basis_->nbf();
  size_t nshells = this->basis_->size();

  //// TODO: Check that P is sane


  auto& tasks = this->load_balancer_->get_tasks();

  //size_t max_npts       = this->load_balancer_->max_npts();
  //size_t max_nbe        = this->load_balancer_->max_nbe();
  //size_t max_npts_x_nbe = this->load_balancer_->max_npts_x_nbe();

  size_t n_deriv = this->func_->is_gga() ? 1 : 0;

  this->timer_.time_op("XCIntegrator.CUDAAlloc", [&](){

    // Allocate Memory
    cuda_data_ = std::make_shared<XCCudaData<value_type>>( );

  });

  // Results
  matrix_type VXC( nbf, nbf );
  value_type  EXC, N_EL;

  this->timer_.time_op("XCIntegrator.LocalWork", [&](){

    // Compute Local contributions to EXC / VXC
    process_batches_cuda_replicated_density_shellbatched_p< value_type>(
      n_deriv, this->timer_, XCWeightAlg::SSF, *this->func_, *this->basis_,
      this->load_balancer_->molecule(), this->load_balancer_->molmeta(),
      *cuda_data_, tasks.begin(), tasks.end(), P.data(), 
      VXC.data(), &EXC, &N_EL 
    );

  } );

  this->timer_.time_op("XCIntegrator.CUDAFree", [&](){
    cuda_data_.reset(); // Free up CUDA memory
  } );

#ifdef GAUXC_ENABLE_MPI


  if( world_size > 1 ) {

    this->timer_.time_op("XCIntegrator.AllReduce", [&]() {

      // Test of communicator is an inter-communicator
      // XXX: Can't think of a case when this would be true, but who knows...
      int inter_flag;
      MPI_Comm_test_inter( this->comm_, &inter_flag );

      // Is Intra-communicator, Allreduce can be done inplace
      if( not inter_flag ) {

        MPI_Allreduce( MPI_IN_PLACE, VXC.data(), nbf*nbf, MPI_DOUBLE,
                       MPI_SUM, this->comm_ );
        MPI_Allreduce( MPI_IN_PLACE, &EXC,  1, MPI_DOUBLE, MPI_SUM, this->comm_ );
        MPI_Allreduce( MPI_IN_PLACE, &N_EL, 1, MPI_DOUBLE, MPI_SUM, this->comm_ );

      // Isn't Intra-communicator (weird), Allreduce can't be done inplace
      } else {

        matrix_type VXC_cpy = VXC;
        value_type EXC_cpy = EXC, N_EL_cpy = N_EL;

        MPI_Allreduce( VXC_cpy.data(), VXC.data(), nbf*nbf, MPI_DOUBLE,
                       MPI_SUM, this->comm_ );
        MPI_Allreduce( &EXC_cpy,  &EXC,  1, MPI_DOUBLE, MPI_SUM, this->comm_ );
        MPI_Allreduce( &N_EL_cpy, &N_EL, 1, MPI_DOUBLE, MPI_SUM, this->comm_ );

      }

    } );

  }

#endif

#ifdef GAUXC_ENABLE_MAGMA
  // Finalize MAGMA
  {
    auto ierr = magma_finalize();
    GAUXC_MAGMA_ERROR( "MAGMA Finalize Failed", ierr );
  }
#endif

  return exc_vxc_type{EXC, std::move(VXC)};

} 

}
}
#endif
