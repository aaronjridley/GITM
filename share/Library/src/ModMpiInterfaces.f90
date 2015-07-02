module ModMPiInterfaces
  implicit none
  private

  public:: mpi_abort
  public:: mpi_allgather
  public:: mpi_allgatherv
  public:: mpi_allreduce
  public:: mpi_barrier
  public:: mpi_bcast
  public:: mpi_bsend
  public:: mpi_comm_create
  public:: mpi_comm_free
  public:: mpi_comm_group
  public:: mpi_comm_rank
  public:: mpi_comm_size
  public:: mpi_file_read
  public:: mpi_file_write
  public:: mpi_finalize
  public:: mpi_gather
  public:: mpi_gatherv
  public:: mpi_group_free
  public:: mpi_group_incl
  public:: mpi_group_range_incl
  public:: mpi_group_rank
  public:: mpi_group_size
  public:: mpi_group_translate_ranks
  public:: mpi_group_union
  public:: mpi_ibsend
  public:: mpi_init
  public:: mpi_irecv
  public:: mpi_irsend
  public:: mpi_isend
  public:: mpi_issend
  public:: mpi_recv
  public:: mpi_reduce
  public:: mpi_rsend
  public:: mpi_send
  public:: mpi_ssend
  public:: mpi_wait
  public:: mpi_waitall


  interface
     subroutine mpi_abort(comm, errorcode, ierror)
       integer, intent(in) :: comm
       integer, intent(in) :: errorcode
       integer, intent(out) :: ierror
     end subroutine mpi_abort
  end interface

  interface mpi_allgather
    module procedure &
    mpi_allgather_i0, &
    mpi_allgather_i1, &
    mpi_allgather_i2, &
    mpi_allgather_r0, &
    mpi_allgather_r1, &
    mpi_allgather_r2, &
    mpi_allgather_r3, &
    mpi_allgather_r4, &
    mpi_allgather_l0, &
    mpi_allgather_l1
  end interface

  interface mpi_allgatherv
    module procedure &
    mpi_allgatherv_i0, &
    mpi_allgatherv_i1, &
    mpi_allgatherv_i2, &
    mpi_allgatherv_r0, &
    mpi_allgatherv_r1, &
    mpi_allgatherv_r2, &
    mpi_allgatherv_r3, &
    mpi_allgatherv_r4, &
    mpi_allgatherv_l0, &
    mpi_allgatherv_l1, &
    mpi_allgatherv_l2
  end interface

  interface mpi_allreduce
    module procedure &
    mpi_allreduce_i0, &
    mpi_allreduce_i1, &
    mpi_allreduce_i2, &
    mpi_allreduce_r0, &
    mpi_allreduce_r1, &
    mpi_allreduce_r2, &
    mpi_allreduce_r3, &
    mpi_allreduce_r4, &
    mpi_allreduce_l0, &
    mpi_allreduce_l1
  end interface

  interface
     subroutine mpi_barrier(comm, ierror) 
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
     end subroutine mpi_barrier
  end interface

  interface mpi_bcast
    module procedure &
    mpi_bcast_i0, &
    mpi_bcast_i1, &
    mpi_bcast_i2, &
    mpi_bcast_r0, &
    mpi_bcast_r1, &
    mpi_bcast_r2, &
    mpi_bcast_r3, &
    mpi_bcast_r4, &
    mpi_bcast_r5, &
    mpi_bcast_s0, &
    mpi_bcast_s1, &
    mpi_bcast_l0, &
    mpi_bcast_l1
  end interface

  interface mpi_bsend
    module procedure &
    mpi_bsend_r0
  end interface

  interface
     subroutine mpi_comm_create(comm, group, newcomm, ierror)
       integer, intent(in) :: comm
       integer, intent(in) :: group
       integer, intent(out) :: newcomm
       integer, intent(out) :: ierror
     end subroutine mpi_comm_create
  end interface

  interface
     subroutine mpi_comm_free(comm, ierror)
       integer, intent(inout) :: comm
       integer, intent(out) :: ierror
     end subroutine mpi_comm_free
  end interface

  interface
     subroutine mpi_comm_group(comm, group, ierror)
       integer, intent(in) :: comm
       integer, intent(out) :: group
       integer, intent(out) :: ierror
     end subroutine mpi_comm_group
  end interface

  interface
     subroutine mpi_comm_rank(comm, rank, ierror)
       integer, intent(in) :: comm
       integer, intent(out) :: rank
       integer, intent(out) :: ierror
     end subroutine mpi_comm_rank
  end interface

  interface
     subroutine mpi_comm_size(comm, size, ierror)
       integer, intent(in) :: comm
       integer, intent(out) :: size
       integer, intent(out) :: ierror
     end subroutine mpi_comm_size
  end interface

  interface mpi_file_read
    module procedure &
    mpi_file_read_i0, &
    mpi_file_read_i1, &
    mpi_file_read_r0, &
    mpi_file_read_r1
  end interface

  interface mpi_file_write
    module procedure &
    mpi_file_write_i0, &
    mpi_file_write_i1, &
    mpi_file_write_r0, &
    mpi_file_write_r1
  end interface

  interface
     subroutine mpi_finalize(ierror)
       integer, intent(out) :: ierror
     end subroutine mpi_finalize
  end interface

  interface mpi_gather
    module procedure &
    mpi_gather_i0, &
    mpi_gather_i1, &
    mpi_gather_i2, &
    mpi_gather_r0, &
    mpi_gather_r1, &
    mpi_gather_r2, &
    mpi_gather_r3, &
    mpi_gather_r4, &
    mpi_gather_l0, &
    mpi_gather_l1, &
    mpi_gather_s0, &
    mpi_gather_s1
  end interface

  interface mpi_gatherv
    module procedure &
    mpi_gatherv_i0, &
    mpi_gatherv_i1, &
    mpi_gatherv_i2, &
    mpi_gatherv_r0, &
    mpi_gatherv_r1, &
    mpi_gatherv_r2, &
    mpi_gatherv_r3, &
    mpi_gatherv_r4, &
    mpi_gatherv_l0, &
    mpi_gatherv_l1
  end interface

  interface
     subroutine mpi_group_free(group, ierror)
       integer, intent(inout) :: group
       integer, intent(out) :: ierror
     end subroutine mpi_group_free
  end interface

  interface
     subroutine mpi_group_incl(group, n, ranks, newgroup, ierror)
       integer, intent(in) :: group
       integer, intent(in) :: n
       integer, intent(in) :: ranks(*)
       integer, intent(out) :: newgroup
       integer, intent(out) :: ierror
     end subroutine mpi_group_incl
  end interface

  interface
     subroutine mpi_group_range_incl(group, n, ranges, newgroup, ierror) 
       integer, intent(in) :: group
       integer, intent(in) :: n
       integer, intent(in) :: ranges(3,*)
       integer, intent(out) :: newgroup
       integer, intent(out) :: ierror
     end subroutine mpi_group_range_incl
  end interface

  interface
     subroutine mpi_group_rank(group, rank, ierror)
       integer, intent(in) :: group
       integer, intent(out) :: rank
       integer, intent(out) :: ierror
     end subroutine mpi_group_rank
  end interface

  interface
     subroutine mpi_group_size(group, size, ierror)
       integer, intent(in) :: group
       integer, intent(out) :: size
       integer, intent(out) :: ierror
     end subroutine mpi_group_size
  end interface

  interface
     subroutine mpi_group_translate_ranks(group1, n, ranks1, &
          group2, ranks2, ierror) 
       integer, intent(in) :: group1
       integer, intent(in) :: n
       integer, intent(in) :: ranks1(*)
       integer, intent(in) :: group2
       integer, intent(out) :: ranks2(*)
       integer, intent(out) :: ierror
     end subroutine mpi_group_translate_ranks
  end interface

  interface
     subroutine mpi_group_union(group1, group2, newgroup, ierror)
       integer, intent(in) :: group1
       integer, intent(in) :: group2
       integer, intent(out) :: newgroup
       integer, intent(out) :: ierror
     end subroutine mpi_group_union
  end interface

  interface mpi_ibsend
    module procedure &
    mpi_ibsend_r0
  end interface

  interface
     subroutine mpi_init(ierror)
       integer, intent(out) :: ierror
     end subroutine mpi_init
  end interface

  interface mpi_irecv
    module procedure &
    mpi_irecv_i0, &
    mpi_irecv_i1, &
    mpi_irecv_i2, &
    mpi_irecv_r0, &
    mpi_irecv_r1, &
    mpi_irecv_r2, &
    mpi_irecv_r3, &
    mpi_irecv_r4, &
    mpi_irecv_l0, &
    mpi_irecv_l1, &
    mpi_irecv_l2
  end interface

  interface mpi_irsend
    module procedure &
    mpi_irsend_r0
  end interface

  interface mpi_isend
    module procedure &
    mpi_isend_i0, &
    mpi_isend_i1, &
    mpi_isend_i2, &
    mpi_isend_r0, &
    mpi_isend_r1, &
    mpi_isend_r2, &
    mpi_isend_r3, &
    mpi_isend_r4, &
    mpi_isend_l0, &
    mpi_isend_l1, &
    mpi_isend_l2
  end interface

  interface mpi_issend
    module procedure &
    mpi_issend_r0
  end interface

  interface mpi_recv
    module procedure &
    mpi_recv_i0, &
    mpi_recv_i1, &
    mpi_recv_i2, &
    mpi_recv_r0, &
    mpi_recv_r1, &
    mpi_recv_r2, &
    mpi_recv_r3, &
    mpi_recv_r4, &
    mpi_recv_l0, &
    mpi_recv_l1, &
    mpi_recv_s0, &
    mpi_recv_s1
  end interface

  interface mpi_reduce
    module procedure &
    mpi_reduce_i0, &
    mpi_reduce_i1, &
    mpi_reduce_i2, &
    mpi_reduce_r0, &
    mpi_reduce_r1, &
    mpi_reduce_r2, &
    mpi_reduce_r3, &
    mpi_reduce_r4, &
    mpi_reduce_l0, &
    mpi_reduce_l1
  end interface

  interface mpi_rsend
    module procedure &
    mpi_rsend_i0, &
    mpi_rsend_i1, &
    mpi_rsend_i2, &
    mpi_rsend_r0, &
    mpi_rsend_r1, &
    mpi_rsend_r2, &
    mpi_rsend_r3, &
    mpi_rsend_r4, &
    mpi_rsend_r5, &
    mpi_rsend_l0, &
    mpi_rsend_l1, &
    mpi_rsend_l2
  end interface

  interface mpi_send
    module procedure &
    mpi_send_i0, &
    mpi_send_i1, &
    mpi_send_i2, &
    mpi_send_r0, &
    mpi_send_r1, &
    mpi_send_r2, &
    mpi_send_r3, &
    mpi_send_r4, &
    mpi_send_s0, &
    mpi_send_s1
  end interface

  interface mpi_ssend
    module procedure &
    mpi_ssend_r0
  end interface

  interface
     subroutine mpi_wait(request, status, ierror)
       use ModMpiOrig, only: mpi_status_size
       integer, intent(inout) :: request
       integer, intent(out) :: status(mpi_status_size)
       integer, intent(out) :: ierror
     end subroutine mpi_wait
  end interface

  interface
     subroutine mpi_waitall(count, array_of_requests, &
          array_of_statuses, ierror) 
       use ModMpiOrig, only: mpi_status_size
       integer, intent(in) :: count
       integer, intent(inout) :: array_of_requests(*)
       integer, intent(out) :: array_of_statuses(mpi_status_size,*)
       integer, intent(out) :: ierror
     end subroutine mpi_waitall
  end interface



contains

     subroutine mpi_allgather_i0(sendbuf, sendcount, sendtype,       &
          recvbuf, recvcount, recvtype, comm, ierror) 
       integer, intent(in) :: sendbuf
       integer, intent(out) :: recvbuf(:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcount
       integer, intent(in) :: recvtype
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allgather

       call mpi_allgather(sendbuf, sendcount, sendtype,       &
          recvbuf, recvcount, recvtype, comm, ierror)
     end subroutine mpi_allgather_i0


     subroutine mpi_allgather_i1(sendbuf, sendcount, sendtype,       &
          recvbuf, recvcount, recvtype, comm, ierror) 
       integer, intent(in) :: sendbuf(:)
       integer, intent(out) :: recvbuf(:,:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcount
       integer, intent(in) :: recvtype
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allgather

       call mpi_allgather(sendbuf, sendcount, sendtype,       &
          recvbuf, recvcount, recvtype, comm, ierror)
     end subroutine mpi_allgather_i1


     subroutine mpi_allgather_i2(sendbuf, sendcount, sendtype,       &
          recvbuf, recvcount, recvtype, comm, ierror) 
       integer, intent(in) :: sendbuf(:,:)
       integer, intent(out) :: recvbuf(:,:,:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcount
       integer, intent(in) :: recvtype
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allgather

       call mpi_allgather(sendbuf, sendcount, sendtype,       &
          recvbuf, recvcount, recvtype, comm, ierror)
     end subroutine mpi_allgather_i2


     subroutine mpi_allgather_r0(sendbuf, sendcount, sendtype,       &
          recvbuf, recvcount, recvtype, comm, ierror) 
       real, intent(in) :: sendbuf
       real, intent(out) :: recvbuf(:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcount
       integer, intent(in) :: recvtype
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allgather

       call mpi_allgather(sendbuf, sendcount, sendtype,       &
          recvbuf, recvcount, recvtype, comm, ierror)
     end subroutine mpi_allgather_r0


     subroutine mpi_allgather_r1(sendbuf, sendcount, sendtype,       &
          recvbuf, recvcount, recvtype, comm, ierror) 
       real, intent(in) :: sendbuf(:)
       real, intent(out) :: recvbuf(:,:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcount
       integer, intent(in) :: recvtype
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allgather

       call mpi_allgather(sendbuf, sendcount, sendtype,       &
          recvbuf, recvcount, recvtype, comm, ierror)
     end subroutine mpi_allgather_r1


     subroutine mpi_allgather_r2(sendbuf, sendcount, sendtype,       &
          recvbuf, recvcount, recvtype, comm, ierror) 
       real, intent(in) :: sendbuf(:,:)
       real, intent(out) :: recvbuf(:,:,:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcount
       integer, intent(in) :: recvtype
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allgather

       call mpi_allgather(sendbuf, sendcount, sendtype,       &
          recvbuf, recvcount, recvtype, comm, ierror)
     end subroutine mpi_allgather_r2


     subroutine mpi_allgather_r3(sendbuf, sendcount, sendtype,       &
          recvbuf, recvcount, recvtype, comm, ierror) 
       real, intent(in) :: sendbuf(:,:,:)
       real, intent(out) :: recvbuf(:,:,:,:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcount
       integer, intent(in) :: recvtype
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allgather

       call mpi_allgather(sendbuf, sendcount, sendtype,       &
          recvbuf, recvcount, recvtype, comm, ierror)
     end subroutine mpi_allgather_r3


     subroutine mpi_allgather_r4(sendbuf, sendcount, sendtype,       &
          recvbuf, recvcount, recvtype, comm, ierror) 
       real, intent(in) :: sendbuf(:,:,:,:)
       real, intent(out) :: recvbuf(:,:,:,:,:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcount
       integer, intent(in) :: recvtype
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allgather

       call mpi_allgather(sendbuf, sendcount, sendtype,       &
          recvbuf, recvcount, recvtype, comm, ierror)
     end subroutine mpi_allgather_r4


     subroutine mpi_allgather_l0(sendbuf, sendcount, sendtype,       &
          recvbuf, recvcount, recvtype, comm, ierror) 
       logical, intent(in) :: sendbuf
       logical, intent(out) :: recvbuf(:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcount
       integer, intent(in) :: recvtype
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allgather

       call mpi_allgather(sendbuf, sendcount, sendtype,       &
          recvbuf, recvcount, recvtype, comm, ierror)
     end subroutine mpi_allgather_l0


     subroutine mpi_allgather_l1(sendbuf, sendcount, sendtype,       &
          recvbuf, recvcount, recvtype, comm, ierror) 
       logical, intent(in) :: sendbuf(:)
       logical, intent(out) :: recvbuf(:,:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcount
       integer, intent(in) :: recvtype
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allgather

       call mpi_allgather(sendbuf, sendcount, sendtype,       &
          recvbuf, recvcount, recvtype, comm, ierror)
     end subroutine mpi_allgather_l1


     subroutine mpi_allgatherv_i0(sendbuf, sendcount, sendtype,      &
          recvbuf, recvcounts, displs, recvtype, comm, ierror) 
       integer, intent(in) :: sendbuf
       integer, intent(out) :: recvbuf
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcounts(:)
       integer, intent(in) :: displs(:)
       integer, intent(in) :: recvtype
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allgatherv

       call mpi_allgatherv(sendbuf, sendcount, sendtype,      &
          recvbuf, recvcounts, displs, recvtype, comm, ierror)
     end subroutine mpi_allgatherv_i0


     subroutine mpi_allgatherv_i1(sendbuf, sendcount, sendtype,      &
          recvbuf, recvcounts, displs, recvtype, comm, ierror) 
       integer, intent(in) :: sendbuf(:)
       integer, intent(out) :: recvbuf(:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcounts(:)
       integer, intent(in) :: displs(:)
       integer, intent(in) :: recvtype
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allgatherv

       call mpi_allgatherv(sendbuf, sendcount, sendtype,      &
          recvbuf, recvcounts, displs, recvtype, comm, ierror)
     end subroutine mpi_allgatherv_i1


     subroutine mpi_allgatherv_i2(sendbuf, sendcount, sendtype,      &
          recvbuf, recvcounts, displs, recvtype, comm, ierror) 
       integer, intent(in) :: sendbuf(:,:)
       integer, intent(out) :: recvbuf(:,:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcounts(:)
       integer, intent(in) :: displs(:)
       integer, intent(in) :: recvtype
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allgatherv

       call mpi_allgatherv(sendbuf, sendcount, sendtype,      &
          recvbuf, recvcounts, displs, recvtype, comm, ierror)
     end subroutine mpi_allgatherv_i2


     subroutine mpi_allgatherv_r0(sendbuf, sendcount, sendtype,      &
          recvbuf, recvcounts, displs, recvtype, comm, ierror) 
       real, intent(in) :: sendbuf
       real, intent(out) :: recvbuf
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcounts(:)
       integer, intent(in) :: displs(:)
       integer, intent(in) :: recvtype
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allgatherv

       call mpi_allgatherv(sendbuf, sendcount, sendtype,      &
          recvbuf, recvcounts, displs, recvtype, comm, ierror)
     end subroutine mpi_allgatherv_r0


     subroutine mpi_allgatherv_r1(sendbuf, sendcount, sendtype,      &
          recvbuf, recvcounts, displs, recvtype, comm, ierror) 
       real, intent(in) :: sendbuf(:)
       real, intent(out) :: recvbuf(:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcounts(:)
       integer, intent(in) :: displs(:)
       integer, intent(in) :: recvtype
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allgatherv

       call mpi_allgatherv(sendbuf, sendcount, sendtype,      &
          recvbuf, recvcounts, displs, recvtype, comm, ierror)
     end subroutine mpi_allgatherv_r1


     subroutine mpi_allgatherv_r2(sendbuf, sendcount, sendtype,      &
          recvbuf, recvcounts, displs, recvtype, comm, ierror) 
       real, intent(in) :: sendbuf(:,:)
       real, intent(out) :: recvbuf(:,:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcounts(:)
       integer, intent(in) :: displs(:)
       integer, intent(in) :: recvtype
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allgatherv

       call mpi_allgatherv(sendbuf, sendcount, sendtype,      &
          recvbuf, recvcounts, displs, recvtype, comm, ierror)
     end subroutine mpi_allgatherv_r2


     subroutine mpi_allgatherv_r3(sendbuf, sendcount, sendtype,      &
          recvbuf, recvcounts, displs, recvtype, comm, ierror) 
       real, intent(in) :: sendbuf(:,:,:)
       real, intent(out) :: recvbuf(:,:,:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcounts(:)
       integer, intent(in) :: displs(:)
       integer, intent(in) :: recvtype
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allgatherv

       call mpi_allgatherv(sendbuf, sendcount, sendtype,      &
          recvbuf, recvcounts, displs, recvtype, comm, ierror)
     end subroutine mpi_allgatherv_r3


     subroutine mpi_allgatherv_r4(sendbuf, sendcount, sendtype,      &
          recvbuf, recvcounts, displs, recvtype, comm, ierror) 
       real, intent(in) :: sendbuf(:,:,:,:)
       real, intent(out) :: recvbuf(:,:,:,:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcounts(:)
       integer, intent(in) :: displs(:)
       integer, intent(in) :: recvtype
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allgatherv

       call mpi_allgatherv(sendbuf, sendcount, sendtype,      &
          recvbuf, recvcounts, displs, recvtype, comm, ierror)
     end subroutine mpi_allgatherv_r4


     subroutine mpi_allgatherv_l0(sendbuf, sendcount, sendtype,      &
          recvbuf, recvcounts, displs, recvtype, comm, ierror) 
       logical, intent(in) :: sendbuf
       logical, intent(out) :: recvbuf
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcounts(:)
       integer, intent(in) :: displs(:)
       integer, intent(in) :: recvtype
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allgatherv

       call mpi_allgatherv(sendbuf, sendcount, sendtype,      &
          recvbuf, recvcounts, displs, recvtype, comm, ierror)
     end subroutine mpi_allgatherv_l0


     subroutine mpi_allgatherv_l1(sendbuf, sendcount, sendtype,      &
          recvbuf, recvcounts, displs, recvtype, comm, ierror) 
       logical, intent(in) :: sendbuf(:)
       logical, intent(out) :: recvbuf(:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcounts(:)
       integer, intent(in) :: displs(:)
       integer, intent(in) :: recvtype
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allgatherv

       call mpi_allgatherv(sendbuf, sendcount, sendtype,      &
          recvbuf, recvcounts, displs, recvtype, comm, ierror)
     end subroutine mpi_allgatherv_l1


     subroutine mpi_allgatherv_l2(sendbuf, sendcount, sendtype,      &
          recvbuf, recvcounts, displs, recvtype, comm, ierror) 
       logical, intent(in) :: sendbuf(:,:)
       logical, intent(out) :: recvbuf(:,:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcounts(:)
       integer, intent(in) :: displs(:)
       integer, intent(in) :: recvtype
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allgatherv

       call mpi_allgatherv(sendbuf, sendcount, sendtype,      &
          recvbuf, recvcounts, displs, recvtype, comm, ierror)
     end subroutine mpi_allgatherv_l2


     subroutine mpi_allreduce_i0(sendbuf, recvbuf, count, datatype,  &
          op, comm, ierror) 
       integer, intent(in) :: sendbuf
       integer, intent(out) :: recvbuf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: op
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allreduce

       call mpi_allreduce(sendbuf, recvbuf, count, datatype,  &
          op, comm, ierror)
     end subroutine mpi_allreduce_i0


     subroutine mpi_allreduce_i1(sendbuf, recvbuf, count, datatype,  &
          op, comm, ierror) 
       integer, intent(in) :: sendbuf(:)
       integer, intent(out) :: recvbuf(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: op
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allreduce

       call mpi_allreduce(sendbuf, recvbuf, count, datatype,  &
          op, comm, ierror)
     end subroutine mpi_allreduce_i1


     subroutine mpi_allreduce_i2(sendbuf, recvbuf, count, datatype,  &
          op, comm, ierror) 
       integer, intent(in) :: sendbuf(:,:)
       integer, intent(out) :: recvbuf(:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: op
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allreduce

       call mpi_allreduce(sendbuf, recvbuf, count, datatype,  &
          op, comm, ierror)
     end subroutine mpi_allreduce_i2


     subroutine mpi_allreduce_r0(sendbuf, recvbuf, count, datatype,  &
          op, comm, ierror) 
       real, intent(in) :: sendbuf
       real, intent(out) :: recvbuf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: op
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allreduce

       call mpi_allreduce(sendbuf, recvbuf, count, datatype,  &
          op, comm, ierror)
     end subroutine mpi_allreduce_r0


     subroutine mpi_allreduce_r1(sendbuf, recvbuf, count, datatype,  &
          op, comm, ierror) 
       real, intent(in) :: sendbuf(:)
       real, intent(out) :: recvbuf(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: op
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allreduce

       call mpi_allreduce(sendbuf, recvbuf, count, datatype,  &
          op, comm, ierror)
     end subroutine mpi_allreduce_r1


     subroutine mpi_allreduce_r2(sendbuf, recvbuf, count, datatype,  &
          op, comm, ierror) 
       real, intent(in) :: sendbuf(:,:)
       real, intent(out) :: recvbuf(:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: op
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allreduce

       call mpi_allreduce(sendbuf, recvbuf, count, datatype,  &
          op, comm, ierror)
     end subroutine mpi_allreduce_r2


     subroutine mpi_allreduce_r3(sendbuf, recvbuf, count, datatype,  &
          op, comm, ierror) 
       real, intent(in) :: sendbuf(:,:,:)
       real, intent(out) :: recvbuf(:,:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: op
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allreduce

       call mpi_allreduce(sendbuf, recvbuf, count, datatype,  &
          op, comm, ierror)
     end subroutine mpi_allreduce_r3


     subroutine mpi_allreduce_r4(sendbuf, recvbuf, count, datatype,  &
          op, comm, ierror) 
       real, intent(in) :: sendbuf(:,:,:,:)
       real, intent(out) :: recvbuf(:,:,:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: op
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allreduce

       call mpi_allreduce(sendbuf, recvbuf, count, datatype,  &
          op, comm, ierror)
     end subroutine mpi_allreduce_r4


     subroutine mpi_allreduce_l0(sendbuf, recvbuf, count, datatype,  &
          op, comm, ierror) 
       logical, intent(in) :: sendbuf
       logical, intent(out) :: recvbuf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: op
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allreduce

       call mpi_allreduce(sendbuf, recvbuf, count, datatype,  &
          op, comm, ierror)
     end subroutine mpi_allreduce_l0


     subroutine mpi_allreduce_l1(sendbuf, recvbuf, count, datatype,  &
          op, comm, ierror) 
       logical, intent(in) :: sendbuf(:)
       logical, intent(out) :: recvbuf(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: op
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_allreduce

       call mpi_allreduce(sendbuf, recvbuf, count, datatype,  &
          op, comm, ierror)
     end subroutine mpi_allreduce_l1


     subroutine mpi_bcast_i0(buffer, count, datatype, root, comm, ierror) 
       integer, intent(inout) :: buffer
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_bcast

       call mpi_bcast(buffer, count, datatype, root, comm, ierror)
     end subroutine mpi_bcast_i0


     subroutine mpi_bcast_i1(buffer, count, datatype, root, comm, ierror) 
       integer, intent(inout) :: buffer(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_bcast

       call mpi_bcast(buffer, count, datatype, root, comm, ierror)
     end subroutine mpi_bcast_i1


     subroutine mpi_bcast_i2(buffer, count, datatype, root, comm, ierror) 
       integer, intent(inout) :: buffer(:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_bcast

       call mpi_bcast(buffer, count, datatype, root, comm, ierror)
     end subroutine mpi_bcast_i2


     subroutine mpi_bcast_r0(buffer, count, datatype, root, comm, ierror) 
       real, intent(inout) :: buffer
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_bcast

       call mpi_bcast(buffer, count, datatype, root, comm, ierror)
     end subroutine mpi_bcast_r0


     subroutine mpi_bcast_r1(buffer, count, datatype, root, comm, ierror) 
       real, intent(inout) :: buffer(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_bcast

       call mpi_bcast(buffer, count, datatype, root, comm, ierror)
     end subroutine mpi_bcast_r1


     subroutine mpi_bcast_r2(buffer, count, datatype, root, comm, ierror) 
       real, intent(inout) :: buffer(:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_bcast

       call mpi_bcast(buffer, count, datatype, root, comm, ierror)
     end subroutine mpi_bcast_r2


     subroutine mpi_bcast_r3(buffer, count, datatype, root, comm, ierror) 
       real, intent(inout) :: buffer(:,:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_bcast

       call mpi_bcast(buffer, count, datatype, root, comm, ierror)
     end subroutine mpi_bcast_r3


     subroutine mpi_bcast_r4(buffer, count, datatype, root, comm, ierror) 
       real, intent(inout) :: buffer(:,:,:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_bcast

       call mpi_bcast(buffer, count, datatype, root, comm, ierror)
     end subroutine mpi_bcast_r4


     subroutine mpi_bcast_r5(buffer, count, datatype, root, comm, ierror) 
       real, intent(inout) :: buffer(:,:,:,:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_bcast

       call mpi_bcast(buffer, count, datatype, root, comm, ierror)
     end subroutine mpi_bcast_r5


     subroutine mpi_bcast_s0(buffer, count, datatype, root, comm, ierror) 
       character(len=*), intent(inout) :: buffer
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_bcast

       call mpi_bcast(buffer, count, datatype, root, comm, ierror)
     end subroutine mpi_bcast_s0


     subroutine mpi_bcast_s1(buffer, count, datatype, root, comm, ierror) 
       character(len=*), intent(inout) :: buffer(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_bcast

       call mpi_bcast(buffer, count, datatype, root, comm, ierror)
     end subroutine mpi_bcast_s1


     subroutine mpi_bcast_l0(buffer, count, datatype, root, comm, ierror) 
       logical, intent(inout) :: buffer
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_bcast

       call mpi_bcast(buffer, count, datatype, root, comm, ierror)
     end subroutine mpi_bcast_l0


     subroutine mpi_bcast_l1(buffer, count, datatype, root, comm, ierror) 
       logical, intent(inout) :: buffer(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_bcast

       call mpi_bcast(buffer, count, datatype, root, comm, ierror)
     end subroutine mpi_bcast_l1


     subroutine mpi_bsend_r0(buf, count, datatype, dest, tag, comm, ierror) 
       real, intent(in) :: buf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_bsend

       call mpi_bsend(buf, count, datatype, dest, tag, comm, ierror)
     end subroutine mpi_bsend_r0


     subroutine mpi_file_read_i0(fh, buf, count, datatype, status, ierror)
       use ModMpiOrig, only: mpi_status_size
       integer, intent(in) :: fh
       integer,  intent(out):: buf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: status(mpi_status_size)
       integer, intent(out) :: ierror
          external mpi_file_read

       call mpi_file_read(fh, buf, count, datatype, status, ierror)
     end subroutine mpi_file_read_i0


     subroutine mpi_file_read_i1(fh, buf, count, datatype, status, ierror)
       use ModMpiOrig, only: mpi_status_size
       integer, intent(in) :: fh
       integer,  intent(out):: buf(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: status(mpi_status_size)
       integer, intent(out) :: ierror
          external mpi_file_read

       call mpi_file_read(fh, buf, count, datatype, status, ierror)
     end subroutine mpi_file_read_i1


     subroutine mpi_file_read_r0(fh, buf, count, datatype, status, ierror)
       use ModMpiOrig, only: mpi_status_size
       integer, intent(in) :: fh
       real,  intent(out):: buf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: status(mpi_status_size)
       integer, intent(out) :: ierror
          external mpi_file_read

       call mpi_file_read(fh, buf, count, datatype, status, ierror)
     end subroutine mpi_file_read_r0


     subroutine mpi_file_read_r1(fh, buf, count, datatype, status, ierror)
       use ModMpiOrig, only: mpi_status_size
       integer, intent(in) :: fh
       real,  intent(out):: buf(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: status(mpi_status_size)
       integer, intent(out) :: ierror
          external mpi_file_read

       call mpi_file_read(fh, buf, count, datatype, status, ierror)
     end subroutine mpi_file_read_r1


     subroutine mpi_file_write_i0(fh, buf, count, datatype, status, ierror)
       use ModMpiOrig, only: mpi_status_size
       integer, intent(in) :: fh
       integer,  intent(in) :: buf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: status(mpi_status_size)
       integer, intent(out) :: ierror
          external mpi_file_write

       call mpi_file_write(fh, buf, count, datatype, status, ierror)
     end subroutine mpi_file_write_i0


     subroutine mpi_file_write_i1(fh, buf, count, datatype, status, ierror)
       use ModMpiOrig, only: mpi_status_size
       integer, intent(in) :: fh
       integer,  intent(in) :: buf(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: status(mpi_status_size)
       integer, intent(out) :: ierror
          external mpi_file_write

       call mpi_file_write(fh, buf, count, datatype, status, ierror)
     end subroutine mpi_file_write_i1


     subroutine mpi_file_write_r0(fh, buf, count, datatype, status, ierror)
       use ModMpiOrig, only: mpi_status_size
       integer, intent(in) :: fh
       real,  intent(in) :: buf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: status(mpi_status_size)
       integer, intent(out) :: ierror
          external mpi_file_write

       call mpi_file_write(fh, buf, count, datatype, status, ierror)
     end subroutine mpi_file_write_r0


     subroutine mpi_file_write_r1(fh, buf, count, datatype, status, ierror)
       use ModMpiOrig, only: mpi_status_size
       integer, intent(in) :: fh
       real,  intent(in) :: buf(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: status(mpi_status_size)
       integer, intent(out) :: ierror
          external mpi_file_write

       call mpi_file_write(fh, buf, count, datatype, status, ierror)
     end subroutine mpi_file_write_r1


     subroutine mpi_gather_i0(sendbuf, sendcount, sendtype, recvbuf, &
          recvcount, recvtype, root, comm, ierror) 
       integer, intent(in) :: sendbuf
       integer, intent(out) :: recvbuf(:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcount
       integer, intent(in) :: recvtype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_gather

       call mpi_gather(sendbuf, sendcount, sendtype, recvbuf, &
          recvcount, recvtype, root, comm, ierror)
     end subroutine mpi_gather_i0


     subroutine mpi_gather_i1(sendbuf, sendcount, sendtype, recvbuf, &
          recvcount, recvtype, root, comm, ierror) 
       integer, intent(in) :: sendbuf(:)
       integer, intent(out) :: recvbuf(:,:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcount
       integer, intent(in) :: recvtype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_gather

       call mpi_gather(sendbuf, sendcount, sendtype, recvbuf, &
          recvcount, recvtype, root, comm, ierror)
     end subroutine mpi_gather_i1


     subroutine mpi_gather_i2(sendbuf, sendcount, sendtype, recvbuf, &
          recvcount, recvtype, root, comm, ierror) 
       integer, intent(in) :: sendbuf(:,:)
       integer, intent(out) :: recvbuf(:,:,:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcount
       integer, intent(in) :: recvtype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_gather

       call mpi_gather(sendbuf, sendcount, sendtype, recvbuf, &
          recvcount, recvtype, root, comm, ierror)
     end subroutine mpi_gather_i2


     subroutine mpi_gather_r0(sendbuf, sendcount, sendtype, recvbuf, &
          recvcount, recvtype, root, comm, ierror) 
       real, intent(in) :: sendbuf
       real, intent(out) :: recvbuf(:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcount
       integer, intent(in) :: recvtype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_gather

       call mpi_gather(sendbuf, sendcount, sendtype, recvbuf, &
          recvcount, recvtype, root, comm, ierror)
     end subroutine mpi_gather_r0


     subroutine mpi_gather_r1(sendbuf, sendcount, sendtype, recvbuf, &
          recvcount, recvtype, root, comm, ierror) 
       real, intent(in) :: sendbuf(:)
       real, intent(out) :: recvbuf(:,:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcount
       integer, intent(in) :: recvtype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_gather

       call mpi_gather(sendbuf, sendcount, sendtype, recvbuf, &
          recvcount, recvtype, root, comm, ierror)
     end subroutine mpi_gather_r1


     subroutine mpi_gather_r2(sendbuf, sendcount, sendtype, recvbuf, &
          recvcount, recvtype, root, comm, ierror) 
       real, intent(in) :: sendbuf(:,:)
       real, intent(out) :: recvbuf(:,:,:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcount
       integer, intent(in) :: recvtype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_gather

       call mpi_gather(sendbuf, sendcount, sendtype, recvbuf, &
          recvcount, recvtype, root, comm, ierror)
     end subroutine mpi_gather_r2


     subroutine mpi_gather_r3(sendbuf, sendcount, sendtype, recvbuf, &
          recvcount, recvtype, root, comm, ierror) 
       real, intent(in) :: sendbuf(:,:,:)
       real, intent(out) :: recvbuf(:,:,:,:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcount
       integer, intent(in) :: recvtype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_gather

       call mpi_gather(sendbuf, sendcount, sendtype, recvbuf, &
          recvcount, recvtype, root, comm, ierror)
     end subroutine mpi_gather_r3


     subroutine mpi_gather_r4(sendbuf, sendcount, sendtype, recvbuf, &
          recvcount, recvtype, root, comm, ierror) 
       real, intent(in) :: sendbuf(:,:,:,:)
       real, intent(out) :: recvbuf(:,:,:,:,:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcount
       integer, intent(in) :: recvtype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_gather

       call mpi_gather(sendbuf, sendcount, sendtype, recvbuf, &
          recvcount, recvtype, root, comm, ierror)
     end subroutine mpi_gather_r4


     subroutine mpi_gather_l0(sendbuf, sendcount, sendtype, recvbuf, &
          recvcount, recvtype, root, comm, ierror) 
       logical, intent(in) :: sendbuf
       logical, intent(out) :: recvbuf(:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcount
       integer, intent(in) :: recvtype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_gather

       call mpi_gather(sendbuf, sendcount, sendtype, recvbuf, &
          recvcount, recvtype, root, comm, ierror)
     end subroutine mpi_gather_l0


     subroutine mpi_gather_l1(sendbuf, sendcount, sendtype, recvbuf, &
          recvcount, recvtype, root, comm, ierror) 
       logical, intent(in) :: sendbuf(:)
       logical, intent(out) :: recvbuf(:,:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcount
       integer, intent(in) :: recvtype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_gather

       call mpi_gather(sendbuf, sendcount, sendtype, recvbuf, &
          recvcount, recvtype, root, comm, ierror)
     end subroutine mpi_gather_l1


     subroutine mpi_gather_s0(sendbuf, sendcount, sendtype, recvbuf, &
          recvcount, recvtype, root, comm, ierror) 
       character(len=*), intent(in) :: sendbuf
       character(len=*), intent(out) :: recvbuf(:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcount
       integer, intent(in) :: recvtype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_gather

       call mpi_gather(sendbuf, sendcount, sendtype, recvbuf, &
          recvcount, recvtype, root, comm, ierror)
     end subroutine mpi_gather_s0


     subroutine mpi_gather_s1(sendbuf, sendcount, sendtype, recvbuf, &
          recvcount, recvtype, root, comm, ierror) 
       character(len=*), intent(in) :: sendbuf(:)
       character(len=*), intent(out) :: recvbuf(:,:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcount
       integer, intent(in) :: recvtype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_gather

       call mpi_gather(sendbuf, sendcount, sendtype, recvbuf, &
          recvcount, recvtype, root, comm, ierror)
     end subroutine mpi_gather_s1


     subroutine mpi_gatherv_i0(sendbuf, sendcount, sendtype,         &
          recvbuf, recvcounts, displs, recvtype, root, comm, ierror) 
       integer, intent(in) :: sendbuf
       integer, intent(out) :: recvbuf
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcounts(:)
       integer, intent(in) :: displs(:)
       integer, intent(in) :: recvtype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_gatherv

       call mpi_gatherv(sendbuf, sendcount, sendtype,         &
          recvbuf, recvcounts, displs, recvtype, root, comm, ierror)
     end subroutine mpi_gatherv_i0


     subroutine mpi_gatherv_i1(sendbuf, sendcount, sendtype,         &
          recvbuf, recvcounts, displs, recvtype, root, comm, ierror) 
       integer, intent(in) :: sendbuf(:)
       integer, intent(out) :: recvbuf(:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcounts(:)
       integer, intent(in) :: displs(:)
       integer, intent(in) :: recvtype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_gatherv

       call mpi_gatherv(sendbuf, sendcount, sendtype,         &
          recvbuf, recvcounts, displs, recvtype, root, comm, ierror)
     end subroutine mpi_gatherv_i1


     subroutine mpi_gatherv_i2(sendbuf, sendcount, sendtype,         &
          recvbuf, recvcounts, displs, recvtype, root, comm, ierror) 
       integer, intent(in) :: sendbuf(:,:)
       integer, intent(out) :: recvbuf(:,:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcounts(:)
       integer, intent(in) :: displs(:)
       integer, intent(in) :: recvtype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_gatherv

       call mpi_gatherv(sendbuf, sendcount, sendtype,         &
          recvbuf, recvcounts, displs, recvtype, root, comm, ierror)
     end subroutine mpi_gatherv_i2


     subroutine mpi_gatherv_r0(sendbuf, sendcount, sendtype,         &
          recvbuf, recvcounts, displs, recvtype, root, comm, ierror) 
       real, intent(in) :: sendbuf
       real, intent(out) :: recvbuf
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcounts(:)
       integer, intent(in) :: displs(:)
       integer, intent(in) :: recvtype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_gatherv

       call mpi_gatherv(sendbuf, sendcount, sendtype,         &
          recvbuf, recvcounts, displs, recvtype, root, comm, ierror)
     end subroutine mpi_gatherv_r0


     subroutine mpi_gatherv_r1(sendbuf, sendcount, sendtype,         &
          recvbuf, recvcounts, displs, recvtype, root, comm, ierror) 
       real, intent(in) :: sendbuf(:)
       real, intent(out) :: recvbuf(:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcounts(:)
       integer, intent(in) :: displs(:)
       integer, intent(in) :: recvtype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_gatherv

       call mpi_gatherv(sendbuf, sendcount, sendtype,         &
          recvbuf, recvcounts, displs, recvtype, root, comm, ierror)
     end subroutine mpi_gatherv_r1


     subroutine mpi_gatherv_r2(sendbuf, sendcount, sendtype,         &
          recvbuf, recvcounts, displs, recvtype, root, comm, ierror) 
       real, intent(in) :: sendbuf(:,:)
       real, intent(out) :: recvbuf(:,:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcounts(:)
       integer, intent(in) :: displs(:)
       integer, intent(in) :: recvtype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_gatherv

       call mpi_gatherv(sendbuf, sendcount, sendtype,         &
          recvbuf, recvcounts, displs, recvtype, root, comm, ierror)
     end subroutine mpi_gatherv_r2


     subroutine mpi_gatherv_r3(sendbuf, sendcount, sendtype,         &
          recvbuf, recvcounts, displs, recvtype, root, comm, ierror) 
       real, intent(in) :: sendbuf(:,:,:)
       real, intent(out) :: recvbuf(:,:,:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcounts(:)
       integer, intent(in) :: displs(:)
       integer, intent(in) :: recvtype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_gatherv

       call mpi_gatherv(sendbuf, sendcount, sendtype,         &
          recvbuf, recvcounts, displs, recvtype, root, comm, ierror)
     end subroutine mpi_gatherv_r3


     subroutine mpi_gatherv_r4(sendbuf, sendcount, sendtype,         &
          recvbuf, recvcounts, displs, recvtype, root, comm, ierror) 
       real, intent(in) :: sendbuf(:,:,:,:)
       real, intent(out) :: recvbuf(:,:,:,:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcounts(:)
       integer, intent(in) :: displs(:)
       integer, intent(in) :: recvtype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_gatherv

       call mpi_gatherv(sendbuf, sendcount, sendtype,         &
          recvbuf, recvcounts, displs, recvtype, root, comm, ierror)
     end subroutine mpi_gatherv_r4


     subroutine mpi_gatherv_l0(sendbuf, sendcount, sendtype,         &
          recvbuf, recvcounts, displs, recvtype, root, comm, ierror) 
       logical, intent(in) :: sendbuf
       logical, intent(out) :: recvbuf
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcounts(:)
       integer, intent(in) :: displs(:)
       integer, intent(in) :: recvtype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_gatherv

       call mpi_gatherv(sendbuf, sendcount, sendtype,         &
          recvbuf, recvcounts, displs, recvtype, root, comm, ierror)
     end subroutine mpi_gatherv_l0


     subroutine mpi_gatherv_l1(sendbuf, sendcount, sendtype,         &
          recvbuf, recvcounts, displs, recvtype, root, comm, ierror) 
       logical, intent(in) :: sendbuf(:)
       logical, intent(out) :: recvbuf(:)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcounts(:)
       integer, intent(in) :: displs(:)
       integer, intent(in) :: recvtype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_gatherv

       call mpi_gatherv(sendbuf, sendcount, sendtype,         &
          recvbuf, recvcounts, displs, recvtype, root, comm, ierror)
     end subroutine mpi_gatherv_l1


     subroutine mpi_ibsend_r0(buf, count, datatype, dest, tag, comm, &
          request, ierror) 
       real, intent(in) :: buf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
          external mpi_ibsend

       call mpi_ibsend(buf, count, datatype, dest, tag, comm, &
          request, ierror)
     end subroutine mpi_ibsend_r0


     subroutine mpi_irecv_i0(buf, count, datatype, source, tag,      &
          comm, request, ierror) 
       integer, intent(out) :: buf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: source
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
          external mpi_irecv

       call mpi_irecv(buf, count, datatype, source, tag,      &
          comm, request, ierror)
     end subroutine mpi_irecv_i0


     subroutine mpi_irecv_i1(buf, count, datatype, source, tag,      &
          comm, request, ierror) 
       integer, intent(out) :: buf(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: source
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
          external mpi_irecv

       call mpi_irecv(buf, count, datatype, source, tag,      &
          comm, request, ierror)
     end subroutine mpi_irecv_i1


     subroutine mpi_irecv_i2(buf, count, datatype, source, tag,      &
          comm, request, ierror) 
       integer, intent(out) :: buf(:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: source
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
          external mpi_irecv

       call mpi_irecv(buf, count, datatype, source, tag,      &
          comm, request, ierror)
     end subroutine mpi_irecv_i2


     subroutine mpi_irecv_r0(buf, count, datatype, source, tag,      &
          comm, request, ierror) 
       real, intent(out) :: buf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: source
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
          external mpi_irecv

       call mpi_irecv(buf, count, datatype, source, tag,      &
          comm, request, ierror)
     end subroutine mpi_irecv_r0


     subroutine mpi_irecv_r1(buf, count, datatype, source, tag,      &
          comm, request, ierror) 
       real, intent(out) :: buf(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: source
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
          external mpi_irecv

       call mpi_irecv(buf, count, datatype, source, tag,      &
          comm, request, ierror)
     end subroutine mpi_irecv_r1


     subroutine mpi_irecv_r2(buf, count, datatype, source, tag,      &
          comm, request, ierror) 
       real, intent(out) :: buf(:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: source
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
          external mpi_irecv

       call mpi_irecv(buf, count, datatype, source, tag,      &
          comm, request, ierror)
     end subroutine mpi_irecv_r2


     subroutine mpi_irecv_r3(buf, count, datatype, source, tag,      &
          comm, request, ierror) 
       real, intent(out) :: buf(:,:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: source
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
          external mpi_irecv

       call mpi_irecv(buf, count, datatype, source, tag,      &
          comm, request, ierror)
     end subroutine mpi_irecv_r3


     subroutine mpi_irecv_r4(buf, count, datatype, source, tag,      &
          comm, request, ierror) 
       real, intent(out) :: buf(:,:,:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: source
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
          external mpi_irecv

       call mpi_irecv(buf, count, datatype, source, tag,      &
          comm, request, ierror)
     end subroutine mpi_irecv_r4


     subroutine mpi_irecv_l0(buf, count, datatype, source, tag,      &
          comm, request, ierror) 
       logical, intent(out) :: buf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: source
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
          external mpi_irecv

       call mpi_irecv(buf, count, datatype, source, tag,      &
          comm, request, ierror)
     end subroutine mpi_irecv_l0


     subroutine mpi_irecv_l1(buf, count, datatype, source, tag,      &
          comm, request, ierror) 
       logical, intent(out) :: buf(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: source
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
          external mpi_irecv

       call mpi_irecv(buf, count, datatype, source, tag,      &
          comm, request, ierror)
     end subroutine mpi_irecv_l1


     subroutine mpi_irecv_l2(buf, count, datatype, source, tag,      &
          comm, request, ierror) 
       logical, intent(out) :: buf(:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: source
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
          external mpi_irecv

       call mpi_irecv(buf, count, datatype, source, tag,      &
          comm, request, ierror)
     end subroutine mpi_irecv_l2


     subroutine mpi_irsend_r0(buf, count, datatype, dest, tag, comm, &
          request, ierror) 
       real, intent(in) :: buf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
          external mpi_irsend

       call mpi_irsend(buf, count, datatype, dest, tag, comm, &
          request, ierror)
     end subroutine mpi_irsend_r0


     subroutine mpi_isend_i0(buf, count, datatype, dest, tag, comm,  &
          request, ierror) 
       integer, intent(in) :: buf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
          external mpi_isend

       call mpi_isend(buf, count, datatype, dest, tag, comm,  &
          request, ierror)
     end subroutine mpi_isend_i0


     subroutine mpi_isend_i1(buf, count, datatype, dest, tag, comm,  &
          request, ierror) 
       integer, intent(in) :: buf(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
          external mpi_isend

       call mpi_isend(buf, count, datatype, dest, tag, comm,  &
          request, ierror)
     end subroutine mpi_isend_i1


     subroutine mpi_isend_i2(buf, count, datatype, dest, tag, comm,  &
          request, ierror) 
       integer, intent(in) :: buf(:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
          external mpi_isend

       call mpi_isend(buf, count, datatype, dest, tag, comm,  &
          request, ierror)
     end subroutine mpi_isend_i2


     subroutine mpi_isend_r0(buf, count, datatype, dest, tag, comm,  &
          request, ierror) 
       real, intent(in) :: buf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
          external mpi_isend

       call mpi_isend(buf, count, datatype, dest, tag, comm,  &
          request, ierror)
     end subroutine mpi_isend_r0


     subroutine mpi_isend_r1(buf, count, datatype, dest, tag, comm,  &
          request, ierror) 
       real, intent(in) :: buf(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
          external mpi_isend

       call mpi_isend(buf, count, datatype, dest, tag, comm,  &
          request, ierror)
     end subroutine mpi_isend_r1


     subroutine mpi_isend_r2(buf, count, datatype, dest, tag, comm,  &
          request, ierror) 
       real, intent(in) :: buf(:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
          external mpi_isend

       call mpi_isend(buf, count, datatype, dest, tag, comm,  &
          request, ierror)
     end subroutine mpi_isend_r2


     subroutine mpi_isend_r3(buf, count, datatype, dest, tag, comm,  &
          request, ierror) 
       real, intent(in) :: buf(:,:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
          external mpi_isend

       call mpi_isend(buf, count, datatype, dest, tag, comm,  &
          request, ierror)
     end subroutine mpi_isend_r3


     subroutine mpi_isend_r4(buf, count, datatype, dest, tag, comm,  &
          request, ierror) 
       real, intent(in) :: buf(:,:,:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
          external mpi_isend

       call mpi_isend(buf, count, datatype, dest, tag, comm,  &
          request, ierror)
     end subroutine mpi_isend_r4


     subroutine mpi_isend_l0(buf, count, datatype, dest, tag, comm,  &
          request, ierror) 
       logical, intent(in) :: buf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
          external mpi_isend

       call mpi_isend(buf, count, datatype, dest, tag, comm,  &
          request, ierror)
     end subroutine mpi_isend_l0


     subroutine mpi_isend_l1(buf, count, datatype, dest, tag, comm,  &
          request, ierror) 
       logical, intent(in) :: buf(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
          external mpi_isend

       call mpi_isend(buf, count, datatype, dest, tag, comm,  &
          request, ierror)
     end subroutine mpi_isend_l1


     subroutine mpi_isend_l2(buf, count, datatype, dest, tag, comm,  &
          request, ierror) 
       logical, intent(in) :: buf(:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
          external mpi_isend

       call mpi_isend(buf, count, datatype, dest, tag, comm,  &
          request, ierror)
     end subroutine mpi_isend_l2


     subroutine mpi_issend_r0(buf, count, datatype, dest, tag, comm, &
          request, ierror) 
       real, intent(in) :: buf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
          external mpi_issend

       call mpi_issend(buf, count, datatype, dest, tag, comm, &
          request, ierror)
     end subroutine mpi_issend_r0


     subroutine mpi_recv_i0(buf, count, datatype, source, tag, comm, &
          status, ierror) 
       use ModMpiOrig, only: mpi_status_size
       integer, intent(out) :: buf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: source
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: status(mpi_status_size)
       integer, intent(out) :: ierror
          external mpi_recv

       call mpi_recv(buf, count, datatype, source, tag, comm, &
          status, ierror)
     end subroutine mpi_recv_i0


     subroutine mpi_recv_i1(buf, count, datatype, source, tag, comm, &
          status, ierror) 
       use ModMpiOrig, only: mpi_status_size
       integer, intent(out) :: buf(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: source
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: status(mpi_status_size)
       integer, intent(out) :: ierror
          external mpi_recv

       call mpi_recv(buf, count, datatype, source, tag, comm, &
          status, ierror)
     end subroutine mpi_recv_i1


     subroutine mpi_recv_i2(buf, count, datatype, source, tag, comm, &
          status, ierror) 
       use ModMpiOrig, only: mpi_status_size
       integer, intent(out) :: buf(:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: source
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: status(mpi_status_size)
       integer, intent(out) :: ierror
          external mpi_recv

       call mpi_recv(buf, count, datatype, source, tag, comm, &
          status, ierror)
     end subroutine mpi_recv_i2


     subroutine mpi_recv_r0(buf, count, datatype, source, tag, comm, &
          status, ierror) 
       use ModMpiOrig, only: mpi_status_size
       real, intent(out) :: buf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: source
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: status(mpi_status_size)
       integer, intent(out) :: ierror
          external mpi_recv

       call mpi_recv(buf, count, datatype, source, tag, comm, &
          status, ierror)
     end subroutine mpi_recv_r0


     subroutine mpi_recv_r1(buf, count, datatype, source, tag, comm, &
          status, ierror) 
       use ModMpiOrig, only: mpi_status_size
       real, intent(out) :: buf(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: source
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: status(mpi_status_size)
       integer, intent(out) :: ierror
          external mpi_recv

       call mpi_recv(buf, count, datatype, source, tag, comm, &
          status, ierror)
     end subroutine mpi_recv_r1


     subroutine mpi_recv_r2(buf, count, datatype, source, tag, comm, &
          status, ierror) 
       use ModMpiOrig, only: mpi_status_size
       real, intent(out) :: buf(:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: source
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: status(mpi_status_size)
       integer, intent(out) :: ierror
          external mpi_recv

       call mpi_recv(buf, count, datatype, source, tag, comm, &
          status, ierror)
     end subroutine mpi_recv_r2


     subroutine mpi_recv_r3(buf, count, datatype, source, tag, comm, &
          status, ierror) 
       use ModMpiOrig, only: mpi_status_size
       real, intent(out) :: buf(:,:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: source
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: status(mpi_status_size)
       integer, intent(out) :: ierror
          external mpi_recv

       call mpi_recv(buf, count, datatype, source, tag, comm, &
          status, ierror)
     end subroutine mpi_recv_r3


     subroutine mpi_recv_r4(buf, count, datatype, source, tag, comm, &
          status, ierror) 
       use ModMpiOrig, only: mpi_status_size
       real, intent(out) :: buf(:,:,:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: source
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: status(mpi_status_size)
       integer, intent(out) :: ierror
          external mpi_recv

       call mpi_recv(buf, count, datatype, source, tag, comm, &
          status, ierror)
     end subroutine mpi_recv_r4


     subroutine mpi_recv_l0(buf, count, datatype, source, tag, comm, &
          status, ierror) 
       use ModMpiOrig, only: mpi_status_size
       logical, intent(out) :: buf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: source
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: status(mpi_status_size)
       integer, intent(out) :: ierror
          external mpi_recv

       call mpi_recv(buf, count, datatype, source, tag, comm, &
          status, ierror)
     end subroutine mpi_recv_l0


     subroutine mpi_recv_l1(buf, count, datatype, source, tag, comm, &
          status, ierror) 
       use ModMpiOrig, only: mpi_status_size
       logical, intent(out) :: buf(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: source
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: status(mpi_status_size)
       integer, intent(out) :: ierror
          external mpi_recv

       call mpi_recv(buf, count, datatype, source, tag, comm, &
          status, ierror)
     end subroutine mpi_recv_l1


     subroutine mpi_recv_s0(buf, count, datatype, source, tag, comm, &
          status, ierror) 
       use ModMpiOrig, only: mpi_status_size
       character(len=*), intent(out) :: buf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: source
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: status(mpi_status_size)
       integer, intent(out) :: ierror
          external mpi_recv

       call mpi_recv(buf, count, datatype, source, tag, comm, &
          status, ierror)
     end subroutine mpi_recv_s0


     subroutine mpi_recv_s1(buf, count, datatype, source, tag, comm, &
          status, ierror) 
       use ModMpiOrig, only: mpi_status_size
       character(len=*), intent(out) :: buf(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: source
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: status(mpi_status_size)
       integer, intent(out) :: ierror
          external mpi_recv

       call mpi_recv(buf, count, datatype, source, tag, comm, &
          status, ierror)
     end subroutine mpi_recv_s1


     subroutine mpi_reduce_i0(sendbuf, recvbuf, count, datatype, op, &
          root, comm, ierror) 
       integer, intent(in) :: sendbuf
       integer, intent(out) :: recvbuf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: op
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_reduce

       call mpi_reduce(sendbuf, recvbuf, count, datatype, op, &
          root, comm, ierror)
     end subroutine mpi_reduce_i0


     subroutine mpi_reduce_i1(sendbuf, recvbuf, count, datatype, op, &
          root, comm, ierror) 
       integer, intent(in) :: sendbuf(:)
       integer, intent(out) :: recvbuf(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: op
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_reduce

       call mpi_reduce(sendbuf, recvbuf, count, datatype, op, &
          root, comm, ierror)
     end subroutine mpi_reduce_i1


     subroutine mpi_reduce_i2(sendbuf, recvbuf, count, datatype, op, &
          root, comm, ierror) 
       integer, intent(in) :: sendbuf(:,:)
       integer, intent(out) :: recvbuf(:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: op
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_reduce

       call mpi_reduce(sendbuf, recvbuf, count, datatype, op, &
          root, comm, ierror)
     end subroutine mpi_reduce_i2


     subroutine mpi_reduce_r0(sendbuf, recvbuf, count, datatype, op, &
          root, comm, ierror) 
       real, intent(in) :: sendbuf
       real, intent(out) :: recvbuf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: op
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_reduce

       call mpi_reduce(sendbuf, recvbuf, count, datatype, op, &
          root, comm, ierror)
     end subroutine mpi_reduce_r0


     subroutine mpi_reduce_r1(sendbuf, recvbuf, count, datatype, op, &
          root, comm, ierror) 
       real, intent(in) :: sendbuf(:)
       real, intent(out) :: recvbuf(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: op
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_reduce

       call mpi_reduce(sendbuf, recvbuf, count, datatype, op, &
          root, comm, ierror)
     end subroutine mpi_reduce_r1


     subroutine mpi_reduce_r2(sendbuf, recvbuf, count, datatype, op, &
          root, comm, ierror) 
       real, intent(in) :: sendbuf(:,:)
       real, intent(out) :: recvbuf(:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: op
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_reduce

       call mpi_reduce(sendbuf, recvbuf, count, datatype, op, &
          root, comm, ierror)
     end subroutine mpi_reduce_r2


     subroutine mpi_reduce_r3(sendbuf, recvbuf, count, datatype, op, &
          root, comm, ierror) 
       real, intent(in) :: sendbuf(:,:,:)
       real, intent(out) :: recvbuf(:,:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: op
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_reduce

       call mpi_reduce(sendbuf, recvbuf, count, datatype, op, &
          root, comm, ierror)
     end subroutine mpi_reduce_r3


     subroutine mpi_reduce_r4(sendbuf, recvbuf, count, datatype, op, &
          root, comm, ierror) 
       real, intent(in) :: sendbuf(:,:,:,:)
       real, intent(out) :: recvbuf(:,:,:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: op
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_reduce

       call mpi_reduce(sendbuf, recvbuf, count, datatype, op, &
          root, comm, ierror)
     end subroutine mpi_reduce_r4


     subroutine mpi_reduce_l0(sendbuf, recvbuf, count, datatype, op, &
          root, comm, ierror) 
       logical, intent(in) :: sendbuf
       logical, intent(out) :: recvbuf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: op
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_reduce

       call mpi_reduce(sendbuf, recvbuf, count, datatype, op, &
          root, comm, ierror)
     end subroutine mpi_reduce_l0


     subroutine mpi_reduce_l1(sendbuf, recvbuf, count, datatype, op, &
          root, comm, ierror) 
       logical, intent(in) :: sendbuf(:)
       logical, intent(out) :: recvbuf(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: op
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_reduce

       call mpi_reduce(sendbuf, recvbuf, count, datatype, op, &
          root, comm, ierror)
     end subroutine mpi_reduce_l1


     subroutine mpi_rsend_i0(buf, count, datatype, dest, tag, comm, ierror) 
       integer, intent(in) :: buf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_rsend

       call mpi_rsend(buf, count, datatype, dest, tag, comm, ierror)
     end subroutine mpi_rsend_i0


     subroutine mpi_rsend_i1(buf, count, datatype, dest, tag, comm, ierror) 
       integer, intent(in) :: buf(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_rsend

       call mpi_rsend(buf, count, datatype, dest, tag, comm, ierror)
     end subroutine mpi_rsend_i1


     subroutine mpi_rsend_i2(buf, count, datatype, dest, tag, comm, ierror) 
       integer, intent(in) :: buf(:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_rsend

       call mpi_rsend(buf, count, datatype, dest, tag, comm, ierror)
     end subroutine mpi_rsend_i2


     subroutine mpi_rsend_r0(buf, count, datatype, dest, tag, comm, ierror) 
       real, intent(in) :: buf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_rsend

       call mpi_rsend(buf, count, datatype, dest, tag, comm, ierror)
     end subroutine mpi_rsend_r0


     subroutine mpi_rsend_r1(buf, count, datatype, dest, tag, comm, ierror) 
       real, intent(in) :: buf(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_rsend

       call mpi_rsend(buf, count, datatype, dest, tag, comm, ierror)
     end subroutine mpi_rsend_r1


     subroutine mpi_rsend_r2(buf, count, datatype, dest, tag, comm, ierror) 
       real, intent(in) :: buf(:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_rsend

       call mpi_rsend(buf, count, datatype, dest, tag, comm, ierror)
     end subroutine mpi_rsend_r2


     subroutine mpi_rsend_r3(buf, count, datatype, dest, tag, comm, ierror) 
       real, intent(in) :: buf(:,:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_rsend

       call mpi_rsend(buf, count, datatype, dest, tag, comm, ierror)
     end subroutine mpi_rsend_r3


     subroutine mpi_rsend_r4(buf, count, datatype, dest, tag, comm, ierror) 
       real, intent(in) :: buf(:,:,:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_rsend

       call mpi_rsend(buf, count, datatype, dest, tag, comm, ierror)
     end subroutine mpi_rsend_r4


     subroutine mpi_rsend_r5(buf, count, datatype, dest, tag, comm, ierror) 
       real, intent(in) :: buf(:,:,:,:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_rsend

       call mpi_rsend(buf, count, datatype, dest, tag, comm, ierror)
     end subroutine mpi_rsend_r5


     subroutine mpi_rsend_l0(buf, count, datatype, dest, tag, comm, ierror) 
       logical, intent(in) :: buf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_rsend

       call mpi_rsend(buf, count, datatype, dest, tag, comm, ierror)
     end subroutine mpi_rsend_l0


     subroutine mpi_rsend_l1(buf, count, datatype, dest, tag, comm, ierror) 
       logical, intent(in) :: buf(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_rsend

       call mpi_rsend(buf, count, datatype, dest, tag, comm, ierror)
     end subroutine mpi_rsend_l1


     subroutine mpi_rsend_l2(buf, count, datatype, dest, tag, comm, ierror) 
       logical, intent(in) :: buf(:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_rsend

       call mpi_rsend(buf, count, datatype, dest, tag, comm, ierror)
     end subroutine mpi_rsend_l2


     subroutine mpi_send_i0(buf, count, datatype, dest, tag, comm, ierror) 
       integer, intent(in) :: buf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_send

       call mpi_send(buf, count, datatype, dest, tag, comm, ierror)
     end subroutine mpi_send_i0


     subroutine mpi_send_i1(buf, count, datatype, dest, tag, comm, ierror) 
       integer, intent(in) :: buf(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_send

       call mpi_send(buf, count, datatype, dest, tag, comm, ierror)
     end subroutine mpi_send_i1


     subroutine mpi_send_i2(buf, count, datatype, dest, tag, comm, ierror) 
       integer, intent(in) :: buf(:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_send

       call mpi_send(buf, count, datatype, dest, tag, comm, ierror)
     end subroutine mpi_send_i2


     subroutine mpi_send_r0(buf, count, datatype, dest, tag, comm, ierror) 
       real, intent(in) :: buf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_send

       call mpi_send(buf, count, datatype, dest, tag, comm, ierror)
     end subroutine mpi_send_r0


     subroutine mpi_send_r1(buf, count, datatype, dest, tag, comm, ierror) 
       real, intent(in) :: buf(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_send

       call mpi_send(buf, count, datatype, dest, tag, comm, ierror)
     end subroutine mpi_send_r1


     subroutine mpi_send_r2(buf, count, datatype, dest, tag, comm, ierror) 
       real, intent(in) :: buf(:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_send

       call mpi_send(buf, count, datatype, dest, tag, comm, ierror)
     end subroutine mpi_send_r2


     subroutine mpi_send_r3(buf, count, datatype, dest, tag, comm, ierror) 
       real, intent(in) :: buf(:,:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_send

       call mpi_send(buf, count, datatype, dest, tag, comm, ierror)
     end subroutine mpi_send_r3


     subroutine mpi_send_r4(buf, count, datatype, dest, tag, comm, ierror) 
       real, intent(in) :: buf(:,:,:,:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_send

       call mpi_send(buf, count, datatype, dest, tag, comm, ierror)
     end subroutine mpi_send_r4


     subroutine mpi_send_s0(buf, count, datatype, dest, tag, comm, ierror) 
       character(len=*), intent(in) :: buf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_send

       call mpi_send(buf, count, datatype, dest, tag, comm, ierror)
     end subroutine mpi_send_s0


     subroutine mpi_send_s1(buf, count, datatype, dest, tag, comm, ierror) 
       character(len=*), intent(in) :: buf(:)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_send

       call mpi_send(buf, count, datatype, dest, tag, comm, ierror)
     end subroutine mpi_send_s1


     subroutine mpi_ssend_r0(buf, count, datatype, dest, tag, comm, ierror) 
       real, intent(in) :: buf
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
          external mpi_ssend

       call mpi_ssend(buf, count, datatype, dest, tag, comm, ierror)
     end subroutine mpi_ssend_r0




end module ModMpiInterfaces
