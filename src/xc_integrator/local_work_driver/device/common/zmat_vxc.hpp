#include "device/xc_device_task.hpp"
#include "device/type_erased_queue.hpp"

namespace GauXC {

void zmat_lda_vxc( size_t        ntasks,
                   int32_t       max_nbf,
                   int32_t       max_npts,
                   XCDeviceTask* tasks_device,
                   type_erased_queue queue );

void zmat_gga_vxc( size_t        ntasks,
                   int32_t       max_nbf,
                   int32_t       max_npts,
                   XCDeviceTask* tasks_device,
                   type_erased_queue queue );

}
