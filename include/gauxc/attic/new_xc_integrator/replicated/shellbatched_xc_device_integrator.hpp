#pragma once
#include <gauxc/gauxc_config.hpp>
#include <gauxc/new_xc_integrator/replicated/replicated_xc_integrator_impl.hpp>
#include <gauxc/new_xc_integrator/xc_integrator_state.hpp>

namespace GauXC  {
namespace detail {

#ifdef GAUXC_ENABLE_DEVICE
template <typename ValueType>
class ShellBatchedXCDeviceIntegrator : public ReplicatedXCIntegratorImpl<ValueType> {

  using base_type  = ReplicatedXCIntegratorImpl<ValueType>;
  using value_type = typename base_type::value_type;
  using basis_type = typename base_type::basis_type;

  XCIntegratorState state_; 

  void eval_exc_vxc_( int64_t m, int64_t n, const value_type* P,
                      int64_t ldp, value_type* VXC, int64_t ldvxc,
                      value_type* EXC ) override;

public:

  template <typename... Args>
  ShellBatchedXCDeviceIntegrator( Args&&... args ) :
    base_type( std::forward<Args>(args)... ) { }

  ShellBatchedXCDeviceIntegrator( const ShellBatchedXCDeviceIntegrator& );
  ShellBatchedXCDeviceIntegrator( ShellBatchedXCDeviceIntegrator&& ) noexcept;

  ~ShellBatchedXCDeviceIntegrator() noexcept;

};

extern template class ShellBatchedXCDeviceIntegrator<double>;
#endif


template <typename ValueType, typename... Args>
std::unique_ptr< ReplicatedXCIntegratorImpl<ValueType> >
  make_shellbatched_device_integrator_impl( Args&&... args ) {

#ifdef GAUXC_ENABLE_DEVICE
  return std::make_unique<ShellBatchedXCDeviceIntegrator<ValueType>>(
    std::forward<Args>(args)...
  );
#else
  std::string msg = std::string(__PRETTY_FUNCTION__)  + 
	            ": GAUXC_ENABLE_DEVICE = FALSE";
  throw std::runtime_error(msg.c_str());
  return nullptr;
#endif

}


}
}
