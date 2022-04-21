#include "xc_device_aos_data.hpp"
#include "buffer_adaptor.hpp"
#include "integrator_util/integrator_common.hpp"
#include <gauxc/exceptions.hpp>

namespace GauXC {

void XCDeviceAoSData::reset_allocations() {
  XCDeviceStackData::reset_allocations(); // Base implementation
  aos_stack.reset();
}

size_t XCDeviceAoSData::get_mem_req( integrator_term_tracker terms,
  const host_task_type& task ) {

  size_t base_size = XCDeviceStackData::get_mem_req(terms, task);
  const auto is_xc_calc = terms.exc_vxc or terms.exc_grad;
  
  // Everything in AoS is not required for current implementations of
  // the weights kernel
  if( not (is_xc_calc or terms.exx) ) return base_size;
  if( is_xc_calc and terms.xc_approx == _UNDEFINED ) 
    GAUXC_GENERIC_EXCEPTION("NO XC APPROX SET FOR XC CALC");

  //const bool is_lda = terms.xc_approx == LDA;
  const bool is_gga = terms.xc_approx == GGA;

  const bool need_grad = is_gga or terms.exc_grad;
  const bool need_hess = is_gga and terms.exc_grad;

  const auto& points           = task.points;
  const auto& submat_cut_bfn   = task.bfn_screening.submat_map;
  const auto& submat_block_bfn = task.bfn_screening.submat_block;
  if( !submat_cut_bfn.size() or !submat_block_bfn.size() )
    GAUXC_GENERIC_EXCEPTION("Must Populate Bfn Submat Maps");

  const auto& submat_cut_cou   = task.cou_screening.submat_map;
  const auto& submat_block_cou = task.cou_screening.submat_block;
  if( terms.exx and (!submat_cut_cou.size() or !submat_block_cou.size()) )
    GAUXC_GENERIC_EXCEPTION("Must Populate Cou Submat Maps");


  // Dimensions
  const size_t npts         = points.size();
  const size_t nbe_bfn      = task.bfn_screening.nbe;
  const size_t ncut_bfn     = submat_cut_bfn.size();
  const size_t nblock_bfn   = submat_block_bfn.size();

  const size_t nbe_cou      = task.cou_screening.nbe;
  const size_t ncut_cou     = submat_cut_cou.size();
  const size_t nblock_cou   = submat_block_cou.size();

  // Collocation + derivatives
  const size_t mem_bf      = nbe_bfn * npts * sizeof(double);
  const size_t mem_bf_grad = need_grad ? 3*mem_bf : 0;
  const size_t mem_bf_hess = need_hess ? 6*mem_bf : 0;

  // LDA/GGA Z Matrix 
  const size_t mem_zmat_lda_gga = is_xc_calc ? nbe_bfn * npts * sizeof(double) : 0;

  // X Matrix graidnet (needed for GGA Gradients)
  const size_t mem_xmat_grad = 
    (is_gga and terms.exc_grad) ? 3 * nbe_bfn * npts * sizeof(double) : 0;

  // F / G Intermediates (needed for EXX calculations)
  const size_t mem_fmat = terms.exx ? nbe_cou * npts * sizeof(double) : 0;
  const size_t mem_gmat = terms.exx ? nbe_cou * npts * sizeof(double) : 0;

  // nbe * nbe scratch
  const size_t mem_nbe_scr = std::max(nbe_cou,nbe_bfn) * nbe_bfn * sizeof(double);

  // Shell index packing (bfn)
  const size_t mem_submat_cut_bfn = 3 * ncut_bfn * sizeof(int32_t);
  const size_t mem_submat_block_bfn = nblock_bfn * sizeof(int32_t);

  // Shell index packing (cou)
  const size_t mem_submat_cut_cou   = terms.exx ? 3 * ncut_cou * sizeof(int32_t) : 0;
  const size_t mem_submat_block_cou = terms.exx ? nblock_cou * sizeof(int32_t)   : 0;



  // Memroty associated with added a task to the indirection
  const size_t mem_task = sizeof(XCDeviceTask);


  return base_size + 
    ( mem_bf + mem_bf_grad + mem_bf_hess + mem_nbe_scr ) +
    ( mem_zmat_lda_gga + mem_xmat_grad )                 +
    ( mem_fmat + mem_gmat )                              +
    ( mem_submat_cut_bfn + mem_submat_block_bfn )        +
    ( mem_submat_cut_cou + mem_submat_block_cou )        +
    ( mem_task );
}





XCDeviceAoSData::device_buffer_t XCDeviceAoSData::allocate_dynamic_stack( 
  integrator_term_tracker terms,
  host_task_iterator task_begin, host_task_iterator task_end, 
  device_buffer_t buf ) {

  // Allocate base info in the stack
  buf = XCDeviceStackData::allocate_dynamic_stack( terms, task_begin, task_end, 
    buf );

  // All data that currently resides in AoS is XC/EXX related and can be skipped
  // for weights
  const auto is_xc_calc = terms.exc_vxc or terms.exc_grad;
  if( not (is_xc_calc or terms.exx) ) return buf; 

  if( is_xc_calc and terms.xc_approx == _UNDEFINED ) 
    GAUXC_GENERIC_EXCEPTION("NO XC APPROX SET FOR XC CALC");
                                               
  const bool is_gga = terms.xc_approx == GGA;

  const bool need_grad = is_gga or  terms.exc_grad;
  const bool need_hess = is_gga and terms.exc_grad;

  // Current Stack
  auto [ ptr, sz ] = buf;
  buffer_adaptor mem( ptr, sz );

  //const size_t submat_chunk_size = this->get_submat_chunk_size(global_dims.nbf,0);


  // Get dimensions
  total_nbe_scr_task_batch   = 0;

  total_nbe_bfn_npts_task_batch = 0; 
  total_ncut_bfn_task_batch     = 0; 
  total_nblock_bfn_task_batch   = 0; 

  total_nbe_cou_npts_task_batch = 0; 
  total_ncut_cou_task_batch     = 0; 
  total_nblock_cou_task_batch   = 0; 
  for( auto it = task_begin; it != task_end; ++it ) {

    const auto& points           = it->points;
    const auto& submat_cut_bfn   = it->bfn_screening.submat_map;
    const auto& submat_block_bfn = it->bfn_screening.submat_block;
    if( !submat_cut_bfn.size() or !submat_block_bfn.size() )
      GAUXC_GENERIC_EXCEPTION("Must Populate Bfn Submat Maps");

    const auto& submat_cut_cou   = it->cou_screening.submat_map;
    const auto& submat_block_cou = it->cou_screening.submat_block;
    if( terms.exx and (!submat_cut_cou.size() or !submat_block_cou.size()) )
      GAUXC_GENERIC_EXCEPTION("Must Populate Cou Submat Maps");

    const size_t npts        = points.size();

    const size_t ncut_bfn    = submat_cut_bfn.size();
    const size_t nblock_bfn  = submat_block_bfn.size();
    const auto nbe_bfn       = it->bfn_screening.nbe;

    const size_t ncut_cou    = submat_cut_cou.size();
    const size_t nblock_cou  = submat_block_cou.size();
    const auto nbe_cou       = it->cou_screening.nbe;

    total_nbe_scr_task_batch      += std::max(nbe_bfn,nbe_cou) * nbe_bfn;

    total_nbe_bfn_npts_task_batch += nbe_bfn * npts;
    total_ncut_bfn_task_batch     += ncut_bfn;
    total_nblock_bfn_task_batch   += nblock_bfn;

    total_nbe_cou_npts_task_batch += nbe_cou * npts;
    total_ncut_cou_task_batch     += ncut_cou;
    total_nblock_cou_task_batch   += nblock_cou;

  }

  // Device task indirection
  const size_t ntask = std::distance( task_begin, task_end );
  aos_stack.device_tasks = mem.aligned_alloc<XCDeviceTask>( ntask , csl);

  // Collocation + derivatives 
  aos_stack.bf_eval_device     = mem.aligned_alloc<double>( total_nbe_bfn_npts_task_batch , csl);
  if( need_grad ) {
    aos_stack.dbf_x_eval_device  = mem.aligned_alloc<double>( total_nbe_bfn_npts_task_batch , csl);
    aos_stack.dbf_y_eval_device  = mem.aligned_alloc<double>( total_nbe_bfn_npts_task_batch , csl);
    aos_stack.dbf_z_eval_device  = mem.aligned_alloc<double>( total_nbe_bfn_npts_task_batch , csl);
  }

  if( need_hess ) {
    aos_stack.d2bf_xx_eval_device  = mem.aligned_alloc<double>( total_nbe_bfn_npts_task_batch , csl);
    aos_stack.d2bf_xy_eval_device  = mem.aligned_alloc<double>( total_nbe_bfn_npts_task_batch , csl);
    aos_stack.d2bf_xz_eval_device  = mem.aligned_alloc<double>( total_nbe_bfn_npts_task_batch , csl);
    aos_stack.d2bf_yy_eval_device  = mem.aligned_alloc<double>( total_nbe_bfn_npts_task_batch , csl);
    aos_stack.d2bf_yz_eval_device  = mem.aligned_alloc<double>( total_nbe_bfn_npts_task_batch , csl);
    aos_stack.d2bf_zz_eval_device  = mem.aligned_alloc<double>( total_nbe_bfn_npts_task_batch , csl);
  }

  // VXC Z Matrix
  if( is_xc_calc ) {
    aos_stack.zmat_vxc_lda_gga_device = mem.aligned_alloc<double>( total_nbe_bfn_npts_task_batch , csl);
  }

  // X Matrix Gradient (for GGA EXC Gradient)
  if( is_gga and terms.exc_grad ) {
    aos_stack.xmat_dx_device = mem.aligned_alloc<double>( total_nbe_bfn_npts_task_batch , csl);
    aos_stack.xmat_dy_device = mem.aligned_alloc<double>( total_nbe_bfn_npts_task_batch , csl);
    aos_stack.xmat_dz_device = mem.aligned_alloc<double>( total_nbe_bfn_npts_task_batch , csl);
  }

  // Scratch buffer
  aos_stack.nbe_scr_device = mem.aligned_alloc<double>( total_nbe_scr_task_batch , csl);

  // AoS buffers (Bfn)
  aos_stack.submat_cut_bfn_device   = mem.aligned_alloc<int32_t>( 3 * total_ncut_bfn_task_batch , csl);
  aos_stack.submat_block_bfn_device = mem.aligned_alloc<int32_t>( total_nblock_bfn_task_batch , csl);

  // EXX Specfic Terms
  if( terms.exx ) {
    aos_stack.fmat_exx_device = mem.aligned_alloc<double>(total_nbe_cou_npts_task_batch, csl);
    aos_stack.gmat_exx_device = mem.aligned_alloc<double>(total_nbe_cou_npts_task_batch, csl);

    // AoS buffers (Bfn)
    aos_stack.submat_cut_cou_device   = mem.aligned_alloc<int32_t>( 3 * total_ncut_cou_task_batch , csl);
    aos_stack.submat_block_cou_device = mem.aligned_alloc<int32_t>( total_nblock_cou_task_batch , csl);
  }


  // Update dynmem data for derived impls
  return device_buffer_t{ mem.stack(), mem.nleft() };
}


void XCDeviceAoSData::pack_and_send( 
  integrator_term_tracker terms,
  host_task_iterator task_begin, host_task_iterator task_end,
  const BasisSetMap& basis_map ) {


  // Pack and send base data
  XCDeviceStackData::pack_and_send( terms, task_begin, task_end, basis_map );

  if( not device_backend_ ) GAUXC_GENERIC_EXCEPTION("Invalid Device Backend");
  const auto is_xc_calc = terms.exc_vxc or terms.exc_grad;

  // All data that currently resides in AoS is XC related and can be skipped
  // for weights
  if( not (is_xc_calc or terms.exx) ) return; 

  if( is_xc_calc and terms.xc_approx == _UNDEFINED ) 
    GAUXC_GENERIC_EXCEPTION("NO XC APPROX SET FOR XC CALC");
  const bool is_lda = terms.xc_approx == LDA;
  const bool is_gga = terms.xc_approx == GGA;

  const bool need_grad = is_gga or terms.exc_grad;
  const bool need_hess = is_gga and terms.exc_grad;

  // Reset AoS
  host_device_tasks.clear();

  // Host Packing Arrays
  std::vector< std::array<int32_t, 3> > submat_cut_bfn_pack;
  std::vector< int32_t > submat_block_bfn_pack;

  std::vector< std::array<int32_t, 3> > submat_cut_cou_pack;
  std::vector< int32_t > submat_block_cou_pack;


  // Contatenation utility
  auto concat_iterable = []( auto& a, const auto& b ) {
    a.insert( a.end(), b.begin(), b.end() );
  };

  //const size_t submat_chunk_size = this->get_submat_chunk_size(global_dims.nbf,0);

  // Pack AoS data and construct indirections
  for( auto it = task_begin; it != task_end; ++it ) {

    const auto  iAtom            = it->iParent;
    const auto& points           = it->points;
    const auto dist_nearest      = it->dist_nearest;

    const auto& submat_cut_bfn   = it->bfn_screening.submat_map;
    const auto& submat_block_bfn = it->bfn_screening.submat_block;
    if( !submat_cut_bfn.size() or !submat_block_bfn.size() )
      GAUXC_GENERIC_EXCEPTION("Must Populate Bfn Submat Maps");

    const auto& submat_cut_cou   = it->cou_screening.submat_map;
    const auto& submat_block_cou = it->cou_screening.submat_block;
    if( terms.exx and (!submat_cut_cou.size() or !submat_block_cou.size()) )
      GAUXC_GENERIC_EXCEPTION("Must Populate Cou Submat Maps");

    // Dimensions
    const size_t npts         = points.size();

    const size_t ncut_bfn     = submat_cut_bfn.size();
    const size_t nblock_bfn   = submat_block_bfn.size();
    const size_t nshells_bfn  = it->bfn_screening.shell_list.size();
    const auto nbe_bfn        = it->bfn_screening.nbe;

    const size_t ncut_cou     = submat_cut_cou.size();
    const size_t nblock_cou   = submat_block_cou.size();
    const size_t nshells_cou  = it->cou_screening.shell_list.size();
    const auto nbe_cou        = it->cou_screening.nbe;


    // Pack Shell indexing
    concat_iterable( submat_cut_bfn_pack, submat_cut_bfn );
    concat_iterable( submat_block_bfn_pack, submat_block_bfn );
    if( terms.exx ) {
      concat_iterable( submat_cut_cou_pack, submat_cut_cou );
      concat_iterable( submat_block_cou_pack, submat_block_cou );
    }

    // Add task to device indirection
    host_device_tasks.emplace_back();

    // Populate indirection with dimensions
    host_device_tasks.back().npts         = npts;
    host_device_tasks.back().iParent      = iAtom;
    host_device_tasks.back().dist_nearest = dist_nearest;

    host_device_tasks.back().bfn_screening.nbe     = nbe_bfn;
    host_device_tasks.back().bfn_screening.ncut    = ncut_bfn;
    host_device_tasks.back().bfn_screening.nblock  = nblock_bfn;
    host_device_tasks.back().bfn_screening.nshells = nshells_bfn;

    if( terms.exx ) {
      host_device_tasks.back().cou_screening.nbe     = nbe_cou;
      host_device_tasks.back().cou_screening.ncut    = ncut_cou;
      host_device_tasks.back().cou_screening.nblock  = nblock_cou;
      host_device_tasks.back().cou_screening.nshells = nshells_cou;
    }

    auto& shell_list_bfn = it->bfn_screening.shell_list;
    host_device_tasks.back().bfn_screening.ibf_begin = 
      basis_map.shell_to_first_ao(shell_list_bfn[0]);

    if( terms.exx ) {
      auto& shell_list_cou = it->cou_screening.shell_list;
      host_device_tasks.back().cou_screening.ibf_begin = 
        basis_map.shell_to_first_ao(shell_list_cou[0]);
    }

  }


  // TODO: Print this if verbose
  //std::cout << "XCDeviceAoSData buf = " << ptr << ", " << sz << std::endl;



  // Send AoS information early to overlap with indirection construction
  device_backend_->copy_async( 3 * submat_cut_bfn_pack.size(), 
    submat_cut_bfn_pack.data()->data(), aos_stack.submat_cut_bfn_device, "send_submat_cut_bfn"  ); 
  device_backend_->copy_async( submat_block_bfn_pack.size(), submat_block_bfn_pack.data(), 
    aos_stack.submat_block_bfn_device, "send_submat_block_bfn"  ); 
  
  if( terms.exx ) {
    device_backend_->copy_async( 3 * submat_cut_cou_pack.size(), 
      submat_cut_cou_pack.data()->data(), aos_stack.submat_cut_cou_device, "send_submat_cut_cou"  ); 
    device_backend_->copy_async( submat_block_cou_pack.size(), submat_block_cou_pack.data(), 
      aos_stack.submat_block_cou_device, "send_submat_block_cou"  ); 
  }

  // Construct full indirection
  {

  const size_t total_npts    = total_npts_task_batch * sizeof(double);
  //buffer_adaptor points_mem ( base_stack.points_device,  3*total_npts );
  buffer_adaptor points_x_mem( base_stack.points_x_device,  total_npts );
  buffer_adaptor points_y_mem( base_stack.points_y_device,  total_npts );
  buffer_adaptor points_z_mem( base_stack.points_z_device,  total_npts );
  buffer_adaptor weights_mem ( base_stack.weights_device,   total_npts );

  const size_t total_ncut_bfn    = total_ncut_bfn_task_batch   * sizeof(int32_t);
  const size_t total_nblock_bfn  = total_nblock_bfn_task_batch * sizeof(int32_t);
  buffer_adaptor submat_cut_bfn_mem( aos_stack.submat_cut_bfn_device, 3*total_ncut_bfn  );
  buffer_adaptor submat_block_bfn_mem( aos_stack.submat_block_bfn_device, total_nblock_bfn);

  const size_t total_ncut_cou    = total_ncut_cou_task_batch   * sizeof(int32_t);
  const size_t total_nblock_cou  = total_nblock_cou_task_batch * sizeof(int32_t);
  buffer_adaptor submat_cut_cou_mem( aos_stack.submat_cut_cou_device, 3*total_ncut_cou  );
  buffer_adaptor submat_block_cou_mem( aos_stack.submat_block_cou_device, total_nblock_cou);

  const size_t total_nbe_scr      = total_nbe_scr_task_batch      * sizeof(double);
  const size_t total_nbe_bfn_npts = total_nbe_bfn_npts_task_batch * sizeof(double);
  const size_t total_nbe_cou_npts = total_nbe_cou_npts_task_batch * sizeof(double);
  buffer_adaptor nbe_mem( aos_stack.nbe_scr_device, total_nbe_scr );
  buffer_adaptor zmat_mem( aos_stack.zmat_vxc_lda_gga_device, total_nbe_bfn_npts );

  buffer_adaptor fmat_mem( aos_stack.fmat_exx_device, total_nbe_cou_npts );
  buffer_adaptor gmat_mem( aos_stack.gmat_exx_device, total_nbe_cou_npts );

  buffer_adaptor bf_mem   ( aos_stack.bf_eval_device,    total_nbe_bfn_npts );
  buffer_adaptor dbf_x_mem( aos_stack.dbf_x_eval_device, total_nbe_bfn_npts );
  buffer_adaptor dbf_y_mem( aos_stack.dbf_y_eval_device, total_nbe_bfn_npts );
  buffer_adaptor dbf_z_mem( aos_stack.dbf_z_eval_device, total_nbe_bfn_npts );

  buffer_adaptor d2bf_xx_mem( aos_stack.d2bf_xx_eval_device, total_nbe_bfn_npts );
  buffer_adaptor d2bf_xy_mem( aos_stack.d2bf_xy_eval_device, total_nbe_bfn_npts );
  buffer_adaptor d2bf_xz_mem( aos_stack.d2bf_xz_eval_device, total_nbe_bfn_npts );
  buffer_adaptor d2bf_yy_mem( aos_stack.d2bf_yy_eval_device, total_nbe_bfn_npts );
  buffer_adaptor d2bf_yz_mem( aos_stack.d2bf_yz_eval_device, total_nbe_bfn_npts );
  buffer_adaptor d2bf_zz_mem( aos_stack.d2bf_zz_eval_device, total_nbe_bfn_npts );

  buffer_adaptor xmat_dx_mem( aos_stack.xmat_dx_device, total_nbe_bfn_npts );
  buffer_adaptor xmat_dy_mem( aos_stack.xmat_dy_device, total_nbe_bfn_npts );
  buffer_adaptor xmat_dz_mem( aos_stack.xmat_dz_device, total_nbe_bfn_npts );

  buffer_adaptor den_mem   ( base_stack.den_eval_device,   total_npts );
  buffer_adaptor dden_x_mem( base_stack.den_x_eval_device, total_npts );
  buffer_adaptor dden_y_mem( base_stack.den_y_eval_device, total_npts );
  buffer_adaptor dden_z_mem( base_stack.den_z_eval_device, total_npts );

  buffer_adaptor eps_mem( base_stack.eps_eval_device, total_npts );
  buffer_adaptor gamma_mem( base_stack.gamma_eval_device, total_npts );
  buffer_adaptor vrho_mem( base_stack.vrho_eval_device, total_npts );
  buffer_adaptor vgamma_mem( base_stack.vgamma_eval_device, total_npts );

  for( auto& task : host_device_tasks ) {
    const auto npts    = task.npts;
    const auto nbe_bfn     = task.bfn_screening.nbe;
    const auto ncut_bfn    = task.bfn_screening.ncut;
    const auto nblock_bfn  = task.bfn_screening.nblock;

    const auto nbe_cou     = task.cou_screening.nbe;
    const auto ncut_cou    = task.cou_screening.ncut;
    const auto nblock_cou  = task.cou_screening.nblock;

    //task.points       = points_mem .aligned_alloc<double>(3*npts, csl);
    task.points_x     = points_x_mem.aligned_alloc<double>(npts, csl);
    task.points_y     = points_y_mem.aligned_alloc<double>(npts, csl);
    task.points_z     = points_z_mem.aligned_alloc<double>(npts, csl);
    task.weights      = weights_mem.aligned_alloc<double>(npts, csl); 

    task.bfn_screening.submat_cut   = submat_cut_bfn_mem.aligned_alloc<int32_t>( 3*ncut_bfn , csl);
    task.bfn_screening.submat_block = submat_block_bfn_mem.aligned_alloc<int32_t>(nblock_bfn, csl);

    task.nbe_scr = nbe_mem .aligned_alloc<double>( std::max(nbe_bfn,nbe_cou) * nbe_bfn, csl);

    if(is_xc_calc) {
      task.zmat = zmat_mem.aligned_alloc<double>( nbe_bfn * npts, csl);
    }

    task.bf = bf_mem.aligned_alloc<double>( nbe_bfn * npts, csl);
    if( need_grad ) {
      task.dbfx = dbf_x_mem.aligned_alloc<double>( nbe_bfn * npts, csl);
      task.dbfy = dbf_y_mem.aligned_alloc<double>( nbe_bfn * npts, csl);
      task.dbfz = dbf_z_mem.aligned_alloc<double>( nbe_bfn * npts, csl);
    }
    if( need_hess ) {
      task.d2bfxx = d2bf_xx_mem.aligned_alloc<double>( nbe_bfn * npts, csl);
      task.d2bfxy = d2bf_xy_mem.aligned_alloc<double>( nbe_bfn * npts, csl);
      task.d2bfxz = d2bf_xz_mem.aligned_alloc<double>( nbe_bfn * npts, csl);
      task.d2bfyy = d2bf_yy_mem.aligned_alloc<double>( nbe_bfn * npts, csl);
      task.d2bfyz = d2bf_yz_mem.aligned_alloc<double>( nbe_bfn * npts, csl);
      task.d2bfzz = d2bf_zz_mem.aligned_alloc<double>( nbe_bfn * npts, csl);
    }

    // X Matrix gradient
    if( is_gga and terms.exc_grad ) {
      task.xmat_x = xmat_dx_mem.aligned_alloc<double>( nbe_bfn * npts, csl);
      task.xmat_y = xmat_dy_mem.aligned_alloc<double>( nbe_bfn * npts, csl);
      task.xmat_z = xmat_dz_mem.aligned_alloc<double>( nbe_bfn * npts, csl);
    }

    if( is_xc_calc ) {
      task.den    = den_mem.aligned_alloc<double>(npts, csl);
      if( need_grad ) {
        task.ddenx  = dden_x_mem.aligned_alloc<double>(npts, csl);
        task.ddeny  = dden_y_mem.aligned_alloc<double>(npts, csl);
        task.ddenz  = dden_z_mem.aligned_alloc<double>(npts, csl);
      }

      task.eps    = eps_mem   .aligned_alloc<double>(npts, csl);
      task.vrho   = vrho_mem  .aligned_alloc<double>(npts, csl);
      if( !is_lda ) {
        task.gamma  = gamma_mem .aligned_alloc<double>(npts, csl);
        task.vgamma = vgamma_mem.aligned_alloc<double>(npts, csl);
      }
    }

    // EXX Specific
    if( terms.exx ) {
      task.cou_screening.submat_cut   = submat_cut_cou_mem.aligned_alloc<int32_t>( 3*ncut_cou , csl);
      task.cou_screening.submat_block = submat_block_cou_mem.aligned_alloc<int32_t>(nblock_cou, csl);
      task.fmat = fmat_mem.aligned_alloc<double>(nbe_cou * npts, csl);
      task.gmat = fmat_mem.aligned_alloc<double>(nbe_cou * npts, csl);
    }

  } // Loop over device tasks

  } // Setup indirection


  // Setup extra pieces to indirection which are algorithm specific
  add_extra_to_indirection(terms, host_device_tasks);

  // Send indirection 
  device_backend_->copy_async( host_device_tasks.size(), host_device_tasks.data(), 
    aos_stack.device_tasks, "send_tasks_device" );


  // Synchronize on the copy stream to keep host vecs in scope
  device_backend_->master_queue_synchronize(); 


}



void XCDeviceAoSData::populate_submat_maps( 
  size_t N,
  host_task_iterator task_begin, host_task_iterator task_end, 
  const BasisSetMap& basis_map ) {


  // Get packing size 
  const size_t submat_chunk_size = this->get_submat_chunk_size(N,0);

  for( auto it = task_begin; it != task_end; ++it ) {

    const auto& shell_list_bfn = it->bfn_screening.shell_list;
    if( shell_list_bfn.size() ) {
      std::tie( it->bfn_screening.submat_map, it->bfn_screening.submat_block ) = 
        gen_compressed_submat_map( basis_map, shell_list_bfn, N, submat_chunk_size );
    }

    const auto& shell_list_cou = it->cou_screening.shell_list;
    if( shell_list_cou.size() ) {
      std::tie( it->cou_screening.submat_map, it->cou_screening.submat_block ) = 
        gen_compressed_submat_map( basis_map, shell_list_cou, N, submat_chunk_size );
    }

  }

}

}
