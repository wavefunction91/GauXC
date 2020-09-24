#include "sycl_zmat.hpp"
#include <gauxc/util/div_ceil.hpp>
#include <gauxc/exceptions/sycl_exception.hpp>

namespace GauXC      {
namespace integrator {
namespace sycl       {

template <typename T>
void zmat_lda_kernel( size_t           ntasks,
                      XCTaskDevice<T>* tasks_device ,
                      cl::sycl::nd_item<3>& item_ct) {

  const size_t batch_idx = item_ct.get_group(2);
  if( batch_idx >= ntasks ) return;

  auto& task = tasks_device[ batch_idx ];
  const auto npts            = task.npts;
  const auto nbf             = task.nbe;
  const auto* vrho_device    = task.vrho;

  const auto* basis_eval_device = task.bf;

  auto* z_matrix_device = task.zmat;

  const size_t tid_x = item_ct.get_global_id(0);
  const size_t tid_y = item_ct.get_global_id(1);

  if( tid_x < nbf and tid_y < npts ) {
    const size_t ibfoff = tid_x + tid_y * nbf;
    const double fact = 0.5 * vrho_device[tid_y];

    z_matrix_device[ tid_x + tid_y * nbf ] = fact * basis_eval_device[ ibfoff ];
  }
}

template <typename T>
void zmat_lda_sycl(size_t ntasks, int32_t max_nbf, int32_t max_npts,
                   XCTaskDevice<T> *tasks_device, cl::sycl::queue *queue) {

  cl::sycl::range<3> threads(16, 16, 1);
  cl::sycl::range<3> blocks(util::div_ceil(max_nbf, threads.get(2)),
                            util::div_ceil(max_npts, threads.get(1)),
                            ntasks);

  GAUXC_SYCL_ERROR( queue->submit([&](cl::sycl::handler &cgh) {
              auto global_range = blocks * threads;

              cgh.parallel_for(cl::sycl::nd_range<3>(global_range, threads),
                               [=](cl::sycl::nd_item<3> item_ct) {
                                   zmat_lda_kernel(ntasks, tasks_device, item_ct);
                               });
          }) );
}

template
void zmat_lda_sycl( size_t                ntasks,
                    int32_t               max_nbf,
                    int32_t               max_npts,
                    XCTaskDevice<double>* tasks_device,
                    cl::sycl::queue*      queue );

template <typename T>
void zmat_gga_kernel( size_t           ntasks,
                      XCTaskDevice<T>* tasks_device ,
                      cl::sycl::nd_item<3>& item_ct) {

  const size_t batch_idx = item_ct.get_group(2);
  if( batch_idx >= ntasks ) return;

  auto& task = tasks_device[ batch_idx ];
  const auto npts            = task.npts;
  const auto nbf             = task.nbe;
  const auto* vrho_device    = task.vrho;
  const auto* vgamma_device  = task.vgamma;
  const auto* den_x_eval_device = task.ddenx;
  const auto* den_y_eval_device = task.ddeny;
  const auto* den_z_eval_device = task.ddenz;

  const auto* basis_eval_device = task.bf;
  const auto* dbasis_x_eval_device = task.dbfx;
  const auto* dbasis_y_eval_device = task.dbfy;
  const auto* dbasis_z_eval_device = task.dbfz;

  auto* z_matrix_device = task.zmat;

  const size_t tid_x = item_ct.get_global_id(0);
  const size_t tid_y = item_ct.get_global_id(1);

  if( tid_x < nbf and tid_y < npts ) {

    const size_t ibfoff = tid_x + tid_y * nbf;
    const double fact_1 = 0.5 * vrho_device[tid_y]  ;
    const double fact_2 = 2.0 * vgamma_device[tid_y];

    const double dx = den_x_eval_device[ tid_y ] * dbasis_x_eval_device[ ibfoff ];
    const double dy = den_y_eval_device[ tid_y ] * dbasis_y_eval_device[ ibfoff ];
    const double dz = den_z_eval_device[ tid_y ] * dbasis_z_eval_device[ ibfoff ];

    z_matrix_device[ tid_x + tid_y * nbf ] =
      fact_1 * basis_eval_device[ ibfoff ] + fact_2 * ( dx + dy + dz );

  }
}

template <typename T>
void zmat_gga_sycl(size_t ntasks, int32_t max_nbf, int32_t max_npts,
                   XCTaskDevice<T> *tasks_device, cl::sycl::queue *queue) {

  cl::sycl::range<3> threads(16, 16, 1);
  cl::sycl::range<3> blocks(util::div_ceil(max_nbf, threads.get(2)),
                            util::div_ceil(max_npts, threads.get(1)),
                            ntasks);

  GAUXC_SYCL_ERROR( queue->submit([&](cl::sycl::handler &cgh) {
          auto global_range = blocks * threads;

          cgh.parallel_for(cl::sycl::nd_range<3>(global_range, threads),
                           [=](cl::sycl::nd_item<3> item_ct) {
                               zmat_gga_kernel(ntasks, tasks_device, item_ct);
                           });
          }) );
}
template
void zmat_gga_sycl( size_t                ntasks,
                    int32_t               max_nbf,
                    int32_t               max_npts,
                    XCTaskDevice<double>* tasks_device,
                    cl::sycl::queue*      queue );

}
}
}
