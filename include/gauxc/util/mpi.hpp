#pragma once
#include <gauxc/gauxc_config.hpp>

#ifdef GAUXC_ENABLE_MPI
#include <mpi.h>
#include <type_traits>

namespace GauXC {

template <typename T>
MPI_Datatype mpi_data_type();

#define REG_MPI_TYPE(TYPE,MPI_TYPE)\
template <> inline MPI_Datatype mpi_data_type<TYPE>(){ return MPI_TYPE; }

REG_MPI_TYPE(int,    MPI_INT     )
REG_MPI_TYPE(double, MPI_DOUBLE  )
REG_MPI_TYPE(size_t, MPI_UINT64_T)

#undef REG_MPI_TYPE

template <typename T>
void allreduce(const T* src, T* dst, int count, MPI_Op op, MPI_Comm comm) {
  MPI_Allreduce(src, dst, count, mpi_data_type<T>(), op, comm);
}


template <typename T>
T allreduce( const T& data, MPI_Op op, MPI_Comm comm) {
  T result;
  allreduce( &data, &result, 1, op, comm);
  return result;
}



template <typename InputIterator, typename OutputIterator>
auto mpi_prefix_sum(InputIterator begin, InputIterator end,
  OutputIterator prefix_sum, MPI_Comm comm) {
  using value_type = typename InputIterator::value_type;
  // Compute local sum
  auto local_sum = std::accumulate(begin, end, value_type(0));

  // Compute global prefix scan (exclusive) to compute local seed values
  value_type prefix_seed = 0;
  MPI_Exscan(&local_sum, &prefix_seed, 1, mpi_data_type<value_type>(),
    MPI_SUM, comm);

  // Compute local exclusive scan
  std::exclusive_scan(begin, end, prefix_sum, value_type(0));

  // Update local scans with seed values
  // XXX: Value on 0 may be clobbered
  int world_rank; MPI_Comm_rank(comm, &world_rank);
  if(world_rank) {
    const size_t n = std::distance(begin,end);
    std::transform(prefix_sum, prefix_sum + n, prefix_sum,
      [=](auto& a){ return a + prefix_seed;});
  }

  return std::make_pair(local_sum, prefix_seed);
}



class MPI_Packed_Buffer {
  MPI_Comm comm_;
  int internal_position_;
  std::vector<std::byte> buffer_;

public:

  MPI_Packed_Buffer(size_t size, MPI_Comm comm) :
    comm_(comm), internal_position_(0), buffer_(size) {}

  auto* buffer() { return buffer_.data(); }
  size_t size()  { return buffer_.size(); }

  template <typename T>
  void pack( const T* ptr, size_t n ) {
    //MPI_Pack( ptr, n, mpi_data_type<T>(), buffer_.data(), buffer_.size(),
    //  &internal_position_, comm_);
    MPI_Pack( ptr, n * sizeof(T), MPI_BYTE, buffer_.data(), buffer_.size(),
      &internal_position_, comm_);
  }

  template <typename T>
  void pack( const T& data ) {
    pack( &data, 1 );
  }


  template <typename T>
  void pack( const std::vector<T>& data ) {
    size_t sz = data.size();
    pack( sz );
    if(sz) pack( data.data(), data.size() );
  }

  template <typename T>
  void unpack( T* ptr, size_t n ) {
    //MPI_Unpack( buffer_.data(), buffer_.size(), &internal_position_,
    //  ptr, n, mpi_data_type<T>(), comm_);
    MPI_Unpack( buffer_.data(), buffer_.size(), &internal_position_,
      ptr, n * sizeof(T), MPI_BYTE, comm_);
  }

  template <typename T>
  void unpack( T& data ) {
    unpack( &data, 1 );
  }

  template <typename T>
  void unpack( std::vector<T>& data ) {
    size_t sz = 0;
    unpack( sz );
    data.resize(sz);
    if(sz) unpack( data.data(), data.size() );
  }

};


template <typename Op>
void ring_execute( const Op& op, MPI_Comm comm ) {
  // Get execution space
  int comm_size, comm_rank;
  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  // Wait for previous rank to send token
  int token;
  if(comm_rank and comm_size > 1)
    MPI_Recv(&token, 1, MPI_INT, comm_rank-1, 0, comm, MPI_STATUS_IGNORE);

  // Execute operation
  op();

  if(comm_size > 1) {
    // Send token to next rank
    MPI_Send(&token, 1, MPI_INT, (comm_rank+1)%comm_size, 0, comm);
    // if Root, wait for final token
    if(!comm_rank)
      MPI_Recv(&token, 1, MPI_INT, comm_size-1, 0, comm, MPI_STATUS_IGNORE);
  }
}


}
#endif
