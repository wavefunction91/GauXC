#include "device/xc_device_task.hpp"

namespace GauXC {

void zmat_lda_vxc_cuda( size_t        ntasks,
                        int32_t       max_nbf,
                        int32_t       max_npts,
                        XCDeviceTask* tasks_device,
                        cudaStream_t  stream );

void zmat_gga_vxc_cuda( size_t        ntasks,
                        int32_t       max_nbf,
                        int32_t       max_npts,
                        XCDeviceTask* tasks_device,
                        cudaStream_t  stream );

}
