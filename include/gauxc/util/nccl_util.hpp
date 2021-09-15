#pragma once
#include <gauxc/gauxc_config.hpp>

#ifdef GAUXC_ENABLE_CUDA
#ifdef GAUXC_ENABLE_MPI
#ifdef GAUXC_ENABLE_NCCL

#include <nccl.h>

namespace GauXC {
namespace util  {

struct nccl_comm {

  ncclComm_t comm;
  
  inline nccl_comm( MPI_Comm mpi_comm ) { 
    int32_t world_rank, world_size;
    MPI_Comm_rank( mpi_comm, &world_rank );
    MPI_Comm_size( mpi_comm, &world_size );

    ncclUniqueId id;
    if (world_rank == 0) ncclGetUniqueId(&id);
    MPI_Bcast(&id, sizeof(id), MPI_BYTE, 0, MPI_COMM_WORLD);

    ncclCommInitRank(&comm, world_size, id, world_rank);
  }

  inline ~nccl_comm() noexcept {
    if( comm != 0 ) ncclCommDestroy(comm);
  }

  nccl_comm( const nccl_comm& ) = delete;
  inline nccl_comm( nccl_comm&& other ) noexcept {
    comm = other.comm;
    other.comm = 0;
  };

  inline operator ncclComm_t() const { return comm; }
};

}
}

#endif
#endif
#endif
