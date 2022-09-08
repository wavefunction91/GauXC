#pragma once
#include "device/scheme1_base.hpp"

#ifdef GAUXC_ENABLE_CUTLASS
#include "cutlass/cutlass.h"
#include "cutlass/gemm/gemm.h"
#endif

namespace GauXC {

struct AoSScheme1CUTLASSBase : public AoSScheme1Base {

  void eval_xmat( XCDeviceData*, bool do_grad ) override final;
  void inc_vxc( XCDeviceData* ) override final;

  struct Data;

  virtual ~AoSScheme1CUTLASSBase() = default;
};

struct AoSScheme1CUTLASSBase::Data : public AoSScheme1Base::Data {

  using base_type = AoSScheme1Base::Data;
  using base_type::host_task_type;
  using base_type::device_buffer_t;

  struct cutlass_data {
    double** dmat_array_device = nullptr;
    double** vmat_array_device = nullptr;
    double** zmat_array_device = nullptr;
    double** bf_array_device   = nullptr;

#ifdef GAUXC_ENABLE_CUTLASS
    cutlass::gemm::GemmCoord* problem_sizes_device = nullptr;
    cutlass::gemm::GemmCoord* syr2k_sizes_device = nullptr;
#endif
    int64_t* ld64_dmat_array_device = nullptr;
    int64_t* ld64_vmat_array_device = nullptr;
    int64_t* ld64_zmat_array_device = nullptr;
    int64_t* ld64_bf_array_device   = nullptr;

    inline void reset(){ std::memset(this,0,sizeof(cutlass_data)); }
  };

  cutlass_data cutlass_stack;

  template <typename... Args>
  Data( Args&&... args ) : base_type( std::forward<Args>(args)... ) { }

  virtual ~Data() = default;

  size_t get_mem_req( integrator_term_tracker, 
    const host_task_type&) override final;
  size_t get_static_mem_requirement() override final; 
  void reset_allocations() override final;
  device_buffer_t allocate_dynamic_stack( integrator_term_tracker terms,
    host_task_iterator begin, host_task_iterator end, device_buffer_t buf )
    override final;
  void pack_and_send( integrator_term_tracker terms,
    host_task_iterator begin, host_task_iterator end, 
    const BasisSetMap& basis_map ) override final;

};

}
