#include <gauxc/util/div_ceil.hpp>

#include "sycl_weights.hpp"
#include "integrator_constants.hpp"
#include "sycl_extensions.hpp"

constexpr double eps_d = std::numeric_limits<double>::epsilon();

namespace GauXC {
    namespace integrator {
        namespace sycl {

            using namespace GauXC::sycl;

            void compute_point_center_dist(size_t npts,
                                           size_t natoms,
                                           const double* coords,
                                           const double* points,
                                           double *dist,
                                           sycl::nd_item<3> item_ct) {

                const int tid_x = item_ct.get_local_id(2) +
                    item_ct.get_group(2) * item_ct.get_local_range().get(2);
                const int tid_y = item_ct.get_local_id(1) +
                    item_ct.get_group(1) * item_ct.get_local_range().get(1);

                if( tid_y < natoms and tid_x < npts ) {
                    const int iAtom = tid_y;
                    const int iPt   = tid_x;

                    const double rx = points[3*iPt + 0] - coords[3*iAtom + 0];
                    const double ry = points[3*iPt + 1] - coords[3*iAtom + 1];
                    const double rz = points[3*iPt + 2] - coords[3*iAtom + 2];

                    dist[iAtom + iPt * natoms] = sycl::sqrt(rx * rx + ry * ry + rz * rz);
                }
            }

            void modify_weights_becke_kernel(size_t npts,
                                             size_t natoms,
                                             const double* RAB,
                                             const double* coords,
                                             const double* dist_scratch,
                                             const int32_t* iparent_device,
                                             double* weights_device,
                                             sycl::nd_item<3> item_ct,
                                             double *shared) {

                auto threadIdx_x = item_ct.get_local_id(2);
                auto threadIdx_y = item_ct.get_local_id(1);

                // Becke partition functions
                auto hBecke = [](double x) {return 1.5 * x - 0.5 * x * x * x;}; // Eq. 19
                auto gBecke = [&](double x) {return hBecke(hBecke(hBecke(x)));}; // Eq. 20 f_3

                for (int ipt = item_ct.get_group(2); ipt < npts;
                     ipt += item_ct.get_group_range(2)) {

                    const auto iParent = iparent_device[ipt];

                    double sum = 0.;
                    double parent_weight = 0.;

                    const double* const local_dist_scratch = dist_scratch + ipt * natoms;
                    for (int iCenter = threadIdx_y; iCenter < natoms;
                         iCenter += item_ct.get_local_range().get(1)) {

                        const double ri = local_dist_scratch[ iCenter ];

                        const double* const local_rab = RAB + iCenter * natoms;

                        double ps = 1.;
                        for (int jCenter = threadIdx_x; jCenter < natoms;
                             jCenter += item_ct.get_local_range().get(2)) {

                            const double rj = local_dist_scratch[ jCenter ];

                            const double mu = (ri - rj) / local_rab[ jCenter ]; // XXX: RAB is symmetric
                            const double s  = 0.5 * ( 1. - gBecke( mu ) );

                            ps *= (iCenter == jCenter) ? 1. : s ;
                        }

                        ps = warp_prod_reduce( item_ct.get_group(), ps ); // XXX: Assumes blockDim.x == 32

                        if( iCenter == iParent ) parent_weight = ps;

                        sum += ps;
                    }

                    // XXX: Assumes blockDim.x == blockDim.y == 32
                    if (threadIdx_x == 0) {
                        shared[threadIdx_y] = sum;
                        shared[threadIdx_y + 1024] = parent_weight;
                    }

                    group_barrier(item_ct.get_group());
                    sum = shared[threadIdx_x];
                    sum = warpReduceSum(sum, item_ct);

                    group_barrier(item_ct.get_group());
                    parent_weight = shared[threadIdx_x + 1024];
                    parent_weight = item_ct.get_sub_group().shuffle(parent_weight, iParent % 32);

                    if (threadIdx_x == 0 and threadIdx_y == 0)
                        weights_device[ipt] *= parent_weight / sum;
                }
            }

            void modify_weights_ssf_kernel(size_t npts,
                                           size_t natoms,
                                           const double* RAB,
                                           const double* coords,
                                           const double* dist_scratch,
                                           const int32_t* iparent_device,
                                           const double* dist_nearest_device,
                                           double* weights_device,
                                           sycl::nd_item<3> item_ct,
                                           double *shared) {

                auto threadIdx_x = item_ct.get_local_id(2);
                auto threadIdx_y = item_ct.get_local_id(1);

                // Frisch partition functions
                auto gFrisch = [](double x) {

                    const double s_x  = x / magic_ssf_factor<>;
                    const double s_x2 = s_x  * s_x;
                    const double s_x3 = s_x  * s_x2;
                    const double s_x5 = s_x3 * s_x2;
                    const double s_x7 = s_x5 * s_x2;

                    return (35.*(s_x - s_x3) + 21.*s_x5 - 5.*s_x7) / 16.;
                };

                auto sFrisch = [&] (double x) {
                    const double g = 0.5 * (1. - gFrisch(x));
                    return (x >= magic_ssf_factor<>) ? 0. : (x <= -magic_ssf_factor<>) ? 1. : g;
                };

                constexpr double weight_tol = 1e-10;

                for (int ipt = item_ct.get_group(2); ipt < npts;
                     ipt += item_ct.get_group_range(2)) {

                    const auto iParent = iparent_device[ipt];

                    double sum = 0.;
                    double parent_weight = 0.;

                    const double* const local_dist_scratch = dist_scratch + ipt * natoms;
                    const double dist_cutoff = 0.5 * (1 - magic_ssf_factor<> ) *
                        dist_nearest_device[ipt];
                    if( local_dist_scratch[iParent] < dist_cutoff ) continue;

                    for (int iCenter = threadIdx_y; iCenter < natoms;
                         iCenter += item_ct.get_local_range().get(1)) {

                        const double ri = local_dist_scratch[ iCenter ];

                        const double* const local_rab = RAB + iCenter * natoms;

                        double ps = 1.;
                        for (int jCenter = threadIdx_x; jCenter < natoms;
                             jCenter += item_ct.get_local_range().get(2))
                            if( sycl::fabs(ps) > weight_tol ) {

                                const double rj = local_dist_scratch[ jCenter ];

                                const double mu = (ri - rj) / local_rab[ jCenter ]; // XXX: RAB is symmetric
                                const double s  = sFrisch( mu );
                                ps *= (iCenter == jCenter) ? 1. : s ;

                            }

                        ps = warp_prod_reduce( item_ct.get_group(), ps ); // XXX: Assumes blockDim.x == 32

                        if( iCenter == iParent ) parent_weight = ps;

                        sum += ps;
                    }

                    // XXX: Assumes blockDim.x == blockDim.y == 32
                    if (threadIdx_x == 0) {
                        shared[threadIdx_y] = sum;
                        shared[threadIdx_y + 1024] = parent_weight;
                    }

                    item_ct.barrier();
                    sum = shared[threadIdx_x];
                    sum = warpReduceSum(sum, item_ct);

                    item_ct.barrier();
                    parent_weight = shared[threadIdx_x + 1024];
                    parent_weight = item_ct.get_sub_group().shuffle(parent_weight, iParent % 32);

                    if (threadIdx_x == 0 and threadIdx_y == 0)
                        weights_device[ipt] *= parent_weight / sum;
                }
            }

// SIMT over points: 1D kernel
            void modify_weights_ssf_kernel_1d(size_t npts,
                                              size_t natoms,
                                              const double* RAB,
                                              const double* coords,
                                              const double* dist_scratch,
                                              const int32_t* iparent_device,
                                              const double* dist_nearest_device,
                                              double* weights_device,
                                              sycl::nd_item<3> item_ct) {

                // Frisch partition functions
                auto gFrisch = [](double x) {

                    const double s_x  = x / magic_ssf_factor<>;
                    const double s_x2 = s_x  * s_x;
                    const double s_x3 = s_x  * s_x2;
                    const double s_x5 = s_x3 * s_x2;
                    const double s_x7 = s_x5 * s_x2;

                    return (35.*(s_x - s_x3) + 21.*s_x5 - 5.*s_x7) / 16.;
                };

#if 0
                auto sFrisch = [&] (double x) {
                    const double g = 0.5 * (1. - gFrisch(x));
                    return (x >= magic_ssf_factor<>) ? 0. : (x <= -magic_ssf_factor<>) ? 1. : g;
                };
#else
                auto sFrisch = [&] (double x) {
                    if( sycl::fabs(x) < magic_ssf_factor<> ) return 0.5 * (1. - gFrisch(x));
                    else if( x >= magic_ssf_factor<> ) return 0.;
                    else                               return 1.;
                };
#endif

                constexpr double weight_tol = 1e-10;

                const int tid_x = item_ct.get_local_id(2) + item_ct.get_group(2) * item_ct.get_local_range().get(2);
                const int nt_x = item_ct.get_local_range().get(2) * item_ct.get_group_range(2);

                //__shared__ double shared[2048];
                for( int ipt = tid_x; ipt < npts; ipt += nt_x ) {

                    const auto iParent = iparent_device[ipt];

                    double sum = 0.;
                    double parent_weight = 0.;

                    const double* const local_dist_scratch = dist_scratch + ipt * natoms;
                    const double dist_cutoff = 0.5 * (1 - magic_ssf_factor<> ) *
                        dist_nearest_device[ipt];
                    if( local_dist_scratch[iParent] < dist_cutoff ) continue;

#if 0
                    for( int iCenter = 0; iCenter < natoms; iCenter++ ) {

                        const double ri = local_dist_scratch[ iCenter ];
                        const double* const local_rab = RAB + iCenter * natoms;

                        double ps = 1.;
                        for( int jCenter = 0; jCenter < natoms; jCenter++ )
                            if( sycl::fabs(ps) > weight_tol ) {
                                if( iCenter != jCenter ) {

                                    const double rj = local_dist_scratch[ jCenter ];
                                    const double mu = (ri - rj) / local_rab[ jCenter ]; // XXX: RAB is symmetric
                                    ps *= sFrisch( mu );

                                }
                            } else break;

                        //__syncwarp();

                        if( iCenter == iParent ) parent_weight = ps;

                        sum += ps;
                    }
#else
                    // Do iParent First
                    {
                        const double ri = local_dist_scratch[ iParent ];
                        const double* const local_rab = RAB + iParent * natoms;

                        parent_weight = 1.;
                        for( int jCenter = 0; jCenter < natoms; jCenter++ )
                            if( parent_weight > weight_tol ) {
                                if( iParent != jCenter ) {

                                    const double rj = local_dist_scratch[ jCenter ];

                                    const double mu = (ri - rj) / local_rab[ jCenter ]; // XXX: RAB is symmetric
                                    parent_weight *= sFrisch( mu );

                                }
                            } else break;

                        //__syncwarp();
                        sum += parent_weight;

                    }

                    if( parent_weight < eps_d ) {
                        weights_device[ipt] = 0.;
                        continue;
                    }

                    for( int iCenter = 0; iCenter < natoms; iCenter++ )
                        if( iParent != iCenter ) {
                            const double ri = local_dist_scratch[ iCenter ];
                            const double* const local_rab = RAB + iCenter * natoms;

                            double ps = 1.;
                            for( int jCenter = 0; jCenter < natoms; jCenter++ )
                                if( ps > weight_tol ) {
                                    if( iCenter != jCenter ) {
                                        const double rj = local_dist_scratch[ jCenter ];

                                        const double mu = (ri - rj) / local_rab[ jCenter ]; // XXX: RAB is symmetric
                                        ps *= sFrisch( mu );
                                    }
                                } else break;

                            //__syncwarp();
                            sum += ps;
                        }
#endif
                    weights_device[ipt] *= parent_weight / sum;
                }
            }

            template <typename F>
            void partition_weights_sycl_SoA(XCWeightAlg weight_alg,
                                            size_t npts, size_t natoms,
                                            const F *points_device,
                                            const int32_t *iparent_device,
                                            const F *dist_nearest_device,
                                            const F *rab_device,
                                            const F *atomic_coords_device,
                                            F *weights_device,
                                            F *dist_scratch_device,
                                            sycl::queue *stream) {

                // Evaluate point-to-atom collocation
                {
                    sycl::range<3> threads(32, 32, 1);
                    sycl::range<3> blocks(util::div_ceil(npts, threads.get(2)),
                                              util::div_ceil(natoms, threads.get(1)));

                    stream->submit([&](sycl::handler &cgh) {
                            auto global_range = blocks * threads;

                            cgh.parallel_for(sycl::nd_range<3>(sycl::range<3>(global_range.get(2),
                                                                                      global_range.get(1),
                                                                                      global_range.get(0)),
                                                                   sycl::range<3>(threads.get(2),
                                                                                      threads.get(1),
                                                                                      threads.get(0))),

                                [=](sycl::nd_item<3> item_ct) {
                                    compute_point_center_dist(npts, natoms, atomic_coords_device,
                                                              points_device, dist_scratch_device,
                                                              item_ct);
                                });
                        });
                }

                const bool partition_weights_1d_kernel = true;

                if( partition_weights_1d_kernel ) {
                    sycl::range<3> threads(1024, 1, 1);
                    sycl::range<3> blocks(util::div_ceil(npts, threads.get(2)));

                    stream->submit([&](sycl::handler &cgh) {
                            auto global_range = blocks * threads;

                            cgh.parallel_for(sycl::nd_range<3>(sycl::range<3>(global_range.get(2),
                                                                                      global_range.get(1),
                                                                                      global_range.get(0)),
                                                                   sycl::range<3>(threads.get(2),
                                                                                      threads.get(1),
                                                                                      threads.get(0))),

                                [=](sycl::nd_item<3> item_ct) {
                                    modify_weights_ssf_kernel_1d(npts, natoms, rab_device, atomic_coords_device,
                                                                 dist_scratch_device, iparent_device, dist_nearest_device,
                                                                 weights_device, item_ct);
                                });
                        });
                }
                else {
                    sycl::range<3> threads(32, 32, 1);
                    sycl::range<3> blocks(npts, 1, 1);
                    auto global_range = blocks * threads;

                    if( weight_alg == XCWeightAlg::SSF )
                    {
                        stream->submit([&](sycl::handler &cgh) {

                                sycl::accessor<double, 1,
                                                   sycl::access::mode::read_write,
                                                   sycl::access::target::local>
                                    shared_acc(sycl::range<1>(2048), cgh);

                                cgh.parallel_for(
                                    sycl::nd_range<3>(sycl::range<3>(global_range.get(2),
                                                                             global_range.get(1),
                                                                             global_range.get(0)),
                                                          sycl::range<3>(threads.get(2),
                                                                             threads.get(1),
                                                                             threads.get(0))),

                                    [=](sycl::nd_item<3> item_ct) {
                                        modify_weights_ssf_kernel(
                                            npts, natoms, rab_device, atomic_coords_device,
                                            dist_scratch_device, iparent_device, dist_nearest_device,
                                            weights_device, item_ct, shared_acc.get_pointer());
                                    });
                            });
                    }
                    else
                    {
                        stream->submit([&](sycl::handler &cgh) {

                                sycl::accessor<double, 1,
                                                   sycl::access::mode::read_write,
                                                   sycl::access::target::local>
                                    shared_acc(sycl::range<1>(2048), cgh);

                                cgh.parallel_for(
                                    sycl::nd_range<3>(sycl::range<3>(global_range.get(2),
                                                                             global_range.get(1),
                                                                             global_range.get(0)),
                                                          sycl::range<3>(threads.get(2),
                                                                             threads.get(1),
                                                                             threads.get(0))),

                                    [=](sycl::nd_item<3> item_ct) {
                                        modify_weights_becke_kernel(
                                            npts, natoms, rab_device, atomic_coords_device,
                                            dist_scratch_device, iparent_device, weights_device, item_ct,
                                            shared_acc.get_pointer());
                                    });
                            });
                    }
                }
            }

            template
            void partition_weights_sycl_SoA( XCWeightAlg    weight_alg,
                                             size_t         npts,
                                             size_t         natoms,
                                             const double*  points_device,
                                             const int32_t* iparent_device,
                                             const double*  dist_nearest_device,
                                             const double*  rab_device,
                                             const double*  atomic_coords_device,
                                             double*  weights_device,
                                             double*  dist_scratch_device,
                                             sycl::queue *stream);

        }
    }
}
