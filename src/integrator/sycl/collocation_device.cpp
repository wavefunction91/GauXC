#include <gauxc/util/div_ceil.hpp>
#include <gauxc/xc_task.hpp>

#include "collocation_petite_kernels.hpp"
#include "collocation_masked_kernels.hpp"
#include "collocation_petite_combined_kernels.hpp"
#include "collocation_masked_combined_kernels.hpp"

namespace GauXC      {
namespace integrator {
namespace sycl       {

    template <typename T>
    void eval_collocation_petite(size_t nshells, size_t nbf, size_t npts,
                                 const Shell<T> *shells_device,
                                 const size_t *offs_device,
                                 const T *pts_device, T *eval_device,
                                 sycl::queue *stream) {

        sycl::range<3> threads(32, 32, 1);
        sycl::range<3> blocks(util::div_ceil(npts, threads[0]),
                                  util::div_ceil(nshells, threads[1]),
                                  1);

        stream->submit([&](sycl::handler &cgh) {
                auto global_range = blocks * threads;

                cgh.parallel_for(
                    sycl::nd_range<3>(sycl::range<3>(global_range.get(2),
                                                             global_range.get(1),
                                                             global_range.get(0)),
                                          sycl::range<3>(threads.get(2),
                                                             threads.get(1),
                                                             threads.get(0))),

                    [=](sycl::nd_item<3> item_ct) {
                        collocation_device_petite_kernel<T>(nshells, nbf, npts, shells_device,
                                                            offs_device, pts_device, eval_device,
                                                            item_ct);
                    });
            });
    }
    template
    void eval_collocation_petite(size_t nshells, size_t nbf, size_t npts,
                                 const Shell<double> *shells_device,
                                 const size_t *offs_device,
                                 const double *pts_device, double *eval_device,
                                 sycl::queue *stream);



    template <typename T>
    void eval_collocation_masked(size_t nshells, size_t nbf, size_t npts,
                                 const Shell<T> *shells_device,
                                 const size_t *mask_device,
                                 const size_t *offs_device,
                                 const T *pts_device, T *eval_device,
                                 sycl::queue *stream) {

        sycl::range<3> threads(32, 32, 1);
        sycl::range<3> blocks(util::div_ceil(npts, threads[0]),
                                  util::div_ceil(nshells, threads[1]), 1);

        stream->submit([&](sycl::handler &cgh) {
                auto global_range = blocks * threads;

                cgh.parallel_for(
                    sycl::nd_range<3>(
                        sycl::range<3>(global_range.get(2), global_range.get(1),
                                           global_range.get(0)),
                        sycl::range<3>(threads.get(2), threads.get(1), threads.get(0))),
                    [=](sycl::nd_item<3> item_ct) {
                        collocation_device_masked_kerne<T>(nshells, nbf, npts, shells_device,
                                                           mask_device, offs_device, pts_device,
                                                           eval_device, item_ct);
                    });
            });
    }
    template
    void eval_collocation_masked(size_t nshells, size_t nbf, size_t npts,
                                 const Shell<double> *shells_device,
                                 const size_t *mask_device,
                                 const size_t *offs_device,
                                 const double *pts_device, double *eval_device,
                                 sycl::queue *stream);



    template <typename T>
    void eval_collocation_petite_combined(size_t ntasks, size_t npts_max,
                                          size_t nshells_max,
                                          XCTaskDevice<T> *device_tasks,
                                          sycl::queue *stream) {

        sycl::range<3> threads(32, 32, 1);
        sycl::range<3> blocks(util::div_ceil(npts_max, threads[0]),
                                  util::div_ceil(nshells_max, threads[1]), ntasks);

        stream->submit([&](sycl::handler &cgh) {
                auto global_range = blocks * threads;

                cgh.parallel_for(
                    sycl::nd_range<3>(
                        sycl::range<3>(global_range.get(2), global_range.get(1),
                                           global_range.get(0)),
                        sycl::range<3>(threads.get(2), threads.get(1), threads.get(0))),
                    [=](sycl::nd_item<3> item_ct) {
                        collocation_device_petite_combined_kernel<T>(ntasks, device_tasks,
                                                                     item_ct);
                    });
            });
    }
    template
    void eval_collocation_petite_combined(size_t ntasks, size_t npts_max,
                                          size_t nshells_max,
                                          XCTaskDevice<double> *device_tasks,
                                          sycl::queue *stream);




    template <typename T>
    void eval_collocation_masked_combined(size_t ntasks, size_t npts_max,
                                          size_t nshells_max,
                                          Shell<T> *shells_device,
                                          XCTaskDevice<T> *device_tasks,
                                          sycl::queue *stream) {

        sycl::range<3> threads(32, 32, 1);
        sycl::range<3> blocks(util::div_ceil(npts_max, threads[0]),
                                  util::div_ceil(nshells_max, threads[1]), ntasks);

        stream->submit([&](sycl::handler &cgh) {
                auto global_range = blocks * threads;

                cgh.parallel_for(
                    sycl::nd_range<3>(
                        sycl::range<3>(global_range.get(2), global_range.get(1),
                                           global_range.get(0)),
                        sycl::range<3>(threads.get(2), threads.get(1), threads.get(0))),
                    [=](sycl::nd_item<3> item_ct) {
                        collocation_device_masked_combined_kernel<T>(ntasks, shells_device,
                                                                     device_tasks, item_ct);
                    });
            });
    }
    template
    void eval_collocation_masked_combined(size_t ntasks, size_t npts_max,
                                          size_t nshells_max,
                                          Shell<double> *shells_device,
                                          XCTaskDevice<double> *device_tasks,
                                          sycl::queue *stream);




    template <typename T>
    void eval_collocation_petite_deriv1(
        size_t nshells, size_t nbf, size_t npts, const Shell<T> *shells_device,
        const size_t *offs_device, const T *pts_device, T *eval_device,
        T *deval_device_x, T *deval_device_y, T *deval_device_z,
        sycl::queue *stream) {

        sycl::range<3> threads(32, 32, 1);
        sycl::range<3> blocks(util::div_ceil(npts, threads[0]),
                                  util::div_ceil(nshells, threads[1]), 1);

        stream->submit([&](sycl::handler &cgh) {
                auto global_range = blocks * threads;

                cgh.parallel_for(
                    sycl::nd_range<3>(
                        sycl::range<3>(global_range.get(2), global_range.get(1),
                                           global_range.get(0)),
                        sycl::range<3>(threads.get(2), threads.get(1), threads.get(0))),
                    [=](sycl::nd_item<3> item_ct) {
                        collocation_device_petite_kernel_deriv1<T>(
                            nshells, nbf, npts, shells_device, offs_device, pts_device,
                            eval_device, deval_device_x, deval_device_y, deval_device_z,
                            item_ct);
                    });
            });
    }
    template
    void eval_collocation_petite_deriv1(
        size_t nshells, size_t nbf, size_t npts, const Shell<double> *shells_device,
        const size_t *offs_device, const double *pts_device, double *eval_device,
        double *deval_device_x, double *deval_device_y, double *deval_device_z,
        sycl::queue *stream);



    template <typename T>
    void eval_collocation_masked_deriv1(
        size_t nshells, size_t nbf, size_t npts, const Shell<T> *shells_device,
        const size_t *mask_device, const size_t *offs_device,
        const T *pts_device, T *eval_device, T *deval_device_x,
        T *deval_device_y, T *deval_device_z, sycl::queue *stream) {

        sycl::range<3> threads(32, 32, 1);
        sycl::range<3> blocks(util::div_ceil(npts, threads[0]),
                                  util::div_ceil(nshells, threads[1]), 1);

        stream->submit([&](sycl::handler &cgh) {
                auto global_range = blocks * threads;

                cgh.parallel_for(
                    sycl::nd_range<3>(
                        sycl::range<3>(global_range.get(2), global_range.get(1),
                                           global_range.get(0)),
                        sycl::range<3>(threads.get(2), threads.get(1), threads.get(0))),
                    [=](sycl::nd_item<3> item_ct) {
                        collocation_device_masked_kernel_deriv1<T>(
                            nshells, nbf, npts, shells_device, mask_device, offs_device,
                            pts_device, eval_device, deval_device_x, deval_device_y,
                            deval_device_z, item_ct);
                    });
            });
    }
    template
    void eval_collocation_masked_deriv1(
        size_t nshells, size_t nbf, size_t npts, const Shell<double> *shells_device,
        const size_t *mask_device, const size_t *offs_device,
        const double *pts_device, double *eval_device, double *deval_device_x,
        double *deval_device_y, double *deval_device_z, sycl::queue *stream);




    template <typename T>
    void eval_collocation_petite_combined_deriv1(size_t ntasks, size_t npts_max,
                                                 size_t nshells_max,
                                                 XCTaskDevice<T> *device_tasks,
                                                 sycl::queue *stream) {

        sycl::range<3> threads(32, 32, 1);
        sycl::range<3> blocks(util::div_ceil(npts_max, threads[0]),
                                  util::div_ceil(nshells_max, threads[1]), ntasks);

        stream->submit([&](sycl::handler &cgh) {
                auto global_range = blocks * threads;

                cgh.parallel_for(
                    sycl::nd_range<3>(
                        sycl::range<3>(global_range.get(2), global_range.get(1),
                                           global_range.get(0)),
                        sycl::range<3>(threads.get(2), threads.get(1), threads.get(0))),
                    [=](sycl::nd_item<3> item_ct) {
                        collocation_device_petite_combined_kernel_deriv1<T>(ntasks, device_tasks,
                                                                            item_ct);
                    });
            });
    }
    template
    void eval_collocation_petite_combined_deriv1(size_t ntasks, size_t npts_max,
                                                 size_t nshells_max,
                                                 XCTaskDevice<double> *device_tasks,
                                                 sycl::queue *stream);



    template <typename T>
    void eval_collocation_masked_combined_deriv1(size_t ntasks, size_t npts_max,
                                                 size_t nshells_max,
                                                 Shell<T> *shells_device,
                                                 XCTaskDevice<T> *device_tasks,
                                                 sycl::queue *stream) {

        sycl::range<3> threads(32, 32, 1);
        sycl::range<3> blocks(util::div_ceil(npts_max, threads[0]),
                                  util::div_ceil(nshells_max, threads[1]),
                                  ntasks);

        stream->submit([&](sycl::handler &cgh) {
                auto global_range = blocks * threads;

                cgh.parallel_for(sycl::nd_range<3>(sycl::range<3>(global_range.get(2),
                                                                          global_range.get(1),
                                                                          global_range.get(0)),
                                                       sycl::range<3>(threads.get(2),
                                                                          threads.get(1),
                                                                          threads.get(0))),

                    [=](sycl::nd_item<3> item_ct) {
                                     collocation_device_masked_combined_kernel_deriv1<T>(
                            ntasks, shells_device, device_tasks, item_ct);
                    });
            });
    }
    template
    void eval_collocation_masked_combined_deriv1(size_t ntasks, size_t npts_max,
                                                 size_t nshells_max,
                                                 Shell<double> *shells_device,
                                                 XCTaskDevice<double> *device_tasks,
                                                 sycl::queue *stream) {


} // namespace sycl
} // namespace integrator
} // namespace GauXC
