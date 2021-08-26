#pragma once

#include <gauxc/new_xc_integrator/impl.hpp>
#include <gauxc/new_xc_integrator/replicated/impl.hpp>
#include <gauxc/new_xc_integrator/replicated/reference_xc_host_integrator.hpp>
#include <gauxc/new_xc_integrator/replicated/incore_xc_device_integrator.hpp>
#include <gauxc/new_xc_integrator/replicated/shellbatched_xc_device_integrator.hpp>

#include <gauxc/util/forward_as_shared_ptr.hpp>

namespace GauXC {

template <typename MatrixType, typename... Args>
XCIntegrator<MatrixType>
  make_default_integrator( ExecutionSpace ex, Args&&... args ) {

  using value_type = typename XCIntegrator<MatrixType>::value_type;

  if( ex == ExecutionSpace::Host ) {

    return XCIntegrator<MatrixType>(
      std::make_unique<detail::ReplicatedXCIntegrator<MatrixType>>(
        detail::make_reference_host_integrator_impl<value_type>( 
          detail::forward_as_shared_ptr(args)... 
        )
      )
    );

  } else {

    return XCIntegrator<MatrixType>(
      std::make_unique<detail::ReplicatedXCIntegrator<MatrixType>>(
        detail::make_incore_device_integrator_impl<value_type>( 
        //detail::make_shellbatched_device_integrator_impl<value_type>( 
          detail::forward_as_shared_ptr(args)... 
        )
      )
    );

  }

}

}
