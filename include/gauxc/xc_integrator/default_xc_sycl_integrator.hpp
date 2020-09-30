#pragma once
#include <gauxc/gauxc_config.hpp>
#include <gauxc/xc_integrator/xc_integrator_impl.hpp>
#include <gauxc/xc_integrator/xc_sycl_data.hpp>
#include <gauxc/xc_integrator/xc_sycl_util.hpp>

#ifdef GAUXC_ENABLE_SYCL
namespace GauXC  {
namespace detail {

using namespace GauXC::integrator::sycl;


template <typename MatrixType>
class DefaultXCSyclIntegrator : public XCIntegratorImpl<MatrixType> {

  using base_type     = XCIntegratorImpl<MatrixType>;
  using matrix_type   = typename base_type::matrix_type;
  using value_type    = typename base_type::value_type;
  using basisset_type = typename base_type::basisset_type;
  using exc_vxc_type  = typename base_type::exc_vxc_type;

  std::shared_ptr< XCSyclData< value_type > > sycl_data_;

  exc_vxc_type eval_exc_vxc_( const MatrixType& ) override;

public:

  template <typename... Args>
  DefaultXCSyclIntegrator( Args&&... args ) :
    base_type( std::forward<Args>(args)... ) { }

  DefaultXCSyclIntegrator( const DefaultXCSyclIntegrator& ) = default;
  DefaultXCSyclIntegrator( DefaultXCSyclIntegrator&& ) noexcept = default;

  ~DefaultXCSyclIntegrator() noexcept = default;

};


template <typename MatrixType>
typename DefaultXCSyclIntegrator<MatrixType>::exc_vxc_type
  DefaultXCSyclIntegrator<MatrixType>::eval_exc_vxc_( const MatrixType& P ) {

  size_t nbf     = this->basis_->nbf();
  size_t nshells = this->basis_->size();

  //// TODO: Check that P is sane


  auto& tasks = this->load_balancer_->get_tasks();

  //size_t max_npts       = this->load_balancer_->max_npts();
  //size_t max_nbe        = this->load_balancer_->max_nbe();
  //size_t max_npts_x_nbe = this->load_balancer_->max_npts_x_nbe();

  size_t n_deriv = this->func_->is_gga() ? 1 : 0;

  // Allocate Memory
  sycl_data_ = std::make_shared<XCSyclData<value_type>>(
    this->load_balancer_->molecule().size(),
    n_deriv,
    nbf,
    nshells,
    false,
    false
  );


  // Results
  matrix_type VXC( nbf, nbf );
  value_type  EXC, N_EL;

  // Compute Local contributions to EXC / VXC
  process_batches_sycl_replicated_p< value_type>(
    n_deriv, XCWeightAlg::SSF, *this->func_, *this->basis_,
    this->load_balancer_->molecule(), this->load_balancer_->molmeta(),
    *sycl_data_, tasks, P.data(), VXC.data(), &EXC, &N_EL
  );

#ifdef GAUXC_ENABLE_MPI

  int world_size;
  MPI_Comm_size( this->comm_, &world_size );

  if( world_size > 1 ) {

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

  }

#endif

  return exc_vxc_type{EXC, std::move(VXC)};

}

}
}
#endif
