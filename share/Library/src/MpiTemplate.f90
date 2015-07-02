!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
# This file is based on the interfaces of the MPI_CHECK library, see
#
# http://andrew.ait.iastate.edu/HPC/MPI-CHECK.htm
#
# The interfaces are slightly modified and simplified. Dimensions for the
# two variable routines are not independent any longer (DIM2 = DIM1 + 1).

module ModMpiTemplate  ! These two lines are here so that 
  interface            ! EMACS can indent the code properly

     subroutine mpi_cancel(request, ierror)
       integer, intent(in) :: request
       integer, intent(out) :: ierror
     end subroutine mpi_cancel

     subroutine mpi_get_count(status, datatype, count, ierror)
       use ModMpiOrig, only: mpi_status_size
       integer, intent(in) :: status(mpi_status_size)
       integer, intent(in) :: datatype
       integer, intent(out) :: count
       integer, intent(out) :: ierror
     end subroutine mpi_get_count

     subroutine mpi_get_elements(status, datatype, count, ierror)
       use ModMpiOrig, only: mpi_status_size
       integer, intent(in) :: status(mpi_status_size)
       integer, intent(in) :: datatype
       integer, intent(out) :: count
       integer, intent(out) :: ierror
     end subroutine mpi_get_elements

     subroutine mpi_iprobe(source, tag, comm, flag, status, ierror)
       use ModMpiOrig, only: mpi_status_size
       logical, intent(out) :: flag
       integer, intent(in) :: source
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: status(mpi_status_size)
       integer, intent(out) :: ierror
     end subroutine mpi_iprobe

     subroutine mpi_pack_size(incount, datatype, comm, size, ierror) 
       integer, intent(in) :: incount
       integer, intent(in) :: datatype
       integer, intent(in) :: comm
       integer, intent(out) :: size
       integer, intent(out) :: ierror
     end subroutine mpi_pack_size

     subroutine mpi_probe(source, tag, comm, status, ierror)
       use ModMpiOrig, only: mpi_status_size
       integer, intent(in) :: source
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: status(mpi_status_size)
       integer, intent(out) :: ierror
     end subroutine mpi_probe

     subroutine mpi_request_free(request, ierror)
       integer, intent(inout) :: request
       integer, intent(out) :: ierror
     end subroutine mpi_request_free

     subroutine mpi_start(request, ierror)
       integer, intent(inout) :: request
       integer, intent(out) :: ierror
     end subroutine mpi_start

     subroutine mpi_startall(count, array_of_requests, ierror)
       integer, intent(in) :: count
       integer, intent(inout) :: array_of_requests(*)
       integer, intent(out) :: ierror
     end subroutine mpi_startall

     subroutine mpi_test(request, flag, status, ierror)
       use ModMpiOrig, only: mpi_status_size
       logical, intent(out) :: flag
       integer, intent(inout) :: request
       integer, intent(out) :: status(mpi_status_size)
       integer, intent(out) :: ierror
     end subroutine mpi_test

     subroutine mpi_testall(count, array_of_requests, flag, &
          array_of_statuses, ierror) 
       use ModMpiOrig, only: mpi_status_size
       logical, intent(out) :: flag
       integer, intent(in) :: count
       integer, intent(inout) :: array_of_requests(*)
       integer, intent(out) :: array_of_statuses(mpi_status_size,*)
       integer, intent(out) :: ierror
     end subroutine mpi_testall

     subroutine mpi_testany(count, array_of_requests, index, flag,  &
          status, ierror) 
       use ModMpiOrig, only: mpi_status_size
       logical flag
       integer, intent(in) :: count
       integer, intent(inout) :: array_of_requests(*)
       integer, intent(out) :: index
       integer status(mpi_status_size)
       integer, intent(out) :: ierror
     end subroutine mpi_testany

     subroutine mpi_testsome(incount, array_of_requests, outcount,  &
          array_of_indices, array_of_statuses, ierror) 
       use ModMpiOrig, only: mpi_status_size
       integer, intent(in) :: incount
       integer, intent(inout) :: array_of_requests(*)
       integer, intent(out) :: outcount
       integer, intent(out) :: array_of_indices(*)
       integer, intent(out) :: array_of_statuses(mpi_status_size,*)
       integer, intent(out) :: ierror
     end subroutine mpi_testsome

     subroutine mpi_test_cancelled(status, flag, ierror)
       use ModMpiOrig, only: mpi_status_size
       logical, intent(out) :: flag
       integer, intent(in) :: status(mpi_status_size)
       integer, intent(out) :: ierror
     end subroutine mpi_test_cancelled

     subroutine mpi_type_commit(datatype, ierror)
       integer, intent(inout) :: datatype
       integer, intent(out) :: ierror
     end subroutine mpi_type_commit

     subroutine mpi_type_contiguous(count, oldtype, newtype, ierror) 
       integer, intent(in) :: count
       integer, intent(in) :: oldtype
       integer, intent(out) :: newtype
       integer, intent(out) :: ierror
     end subroutine mpi_type_contiguous

     subroutine mpi_type_extent(datatype, extent, ierror)
       integer, intent(in) :: datatype
       integer, intent(out) :: extent
       integer, intent(out) :: ierror
     end subroutine mpi_type_extent

     subroutine mpi_type_free(datatype, ierror)
       integer, intent(inout) :: datatype
       integer, intent(out) :: ierror
     end subroutine mpi_type_free

     subroutine mpi_type_hindexed(count, array_of_blocklengths, &
          array_of_displacements, oldtype, newtype, ierror) 
       integer, intent(in) :: count
       integer, intent(in) :: array_of_blocklengths(*)
       integer, intent(in) :: array_of_displacements(*)
       integer, intent(in) :: oldtype
       integer, intent(out) :: newtype
       integer, intent(out) :: ierror
     end subroutine mpi_type_hindexed

     subroutine mpi_type_hvector(count, blocklength, stride, &
          oldtype, newtype, ierror) 
       integer, intent(in) :: count
       integer, intent(in) :: blocklength
       integer, intent(in) :: stride
       integer, intent(in) :: oldtype
       integer, intent(out) :: newtype
       integer, intent(out) :: ierror
     end subroutine mpi_type_hvector

     subroutine mpi_type_indexed(count, array_of_blocklengths, &
          array_of_displacements, oldtype, newtype, ierror) 
       integer, intent(in) :: count
       integer, intent(in) :: array_of_blocklengths(*)
       integer, intent(in) :: array_of_displacements(*)
       integer, intent(in) :: oldtype
       integer, intent(out) :: newtype
       integer, intent(out) :: ierror
     end subroutine mpi_type_indexed

     subroutine mpi_type_lb( datatype, displacement, ierror)
       integer, intent(in) :: datatype
       integer, intent(out) :: displacement
       integer, intent(out) :: ierror
     end subroutine mpi_type_lb

     subroutine mpi_type_size(datatype, size, ierror)
       integer, intent(in) :: datatype
       integer, intent(out) :: size
       integer, intent(out) :: ierror
     end subroutine mpi_type_size

     subroutine mpi_type_struct(count, array_of_blocklengths, &
          array_of_displacements, array_of_types, newtype, ierror) 
       integer, intent(in) :: count
       integer array_of_blocklengths(*)
       integer, intent(in) :: array_of_displacements(*)
       integer, intent(in) :: array_of_types(*)
       integer, intent(out) :: newtype
       integer, intent(out) :: ierror
     end subroutine mpi_type_struct

     subroutine mpi_type_ub( datatype, displacement, ierror)
       integer, intent(in) :: datatype
       integer, intent(out) :: displacement
       integer, intent(out) :: ierror
     end subroutine mpi_type_ub

     subroutine mpi_type_vector(count, blocklength, stride, &
          oldtype, newtype, ierror) 
       integer, intent(in) :: count
       integer, intent(in) :: blocklength
       integer, intent(in) :: stride
       integer, intent(in) :: oldtype
       integer, intent(out) :: newtype
       integer, intent(out) :: ierror
     end subroutine mpi_type_vector

     subroutine mpi_wait(request, status, ierror)
       use ModMpiOrig, only: mpi_status_size
       integer, intent(inout) :: request
       integer, intent(out) :: status(mpi_status_size)
       integer, intent(out) :: ierror
     end subroutine mpi_wait

     subroutine mpi_waitall(count, array_of_requests, &
          array_of_statuses, ierror) 
       use ModMpiOrig, only: mpi_status_size
       integer, intent(in) :: count
       integer, intent(inout) :: array_of_requests(*)
       integer, intent(out) :: array_of_statuses(mpi_status_size,*)
       integer, intent(out) :: ierror
     end subroutine mpi_waitall

     subroutine mpi_waitany(count, array_of_requests, index, status, ierror) 
       use ModMpiOrig, only: mpi_status_size
       integer, intent(in) :: count
       integer, intent(inout) :: array_of_requests(*)
       integer, intent(out) :: index
       integer, intent(out) :: status(mpi_status_size)
       integer, intent(out) :: ierror
     end subroutine mpi_waitany

     subroutine mpi_waitsome(incount, array_of_requests, outcount,  &
          array_of_indices, array_of_statuses, ierror) 
       use ModMpiOrig, only: mpi_status_size
       integer, intent(in) :: incount
       integer, intent(inout) :: array_of_requests(*)
       integer, intent(out) :: outcount
       integer, intent(out) :: array_of_indices(*)
       integer, intent(out) :: array_of_statuses(mpi_status_size,*)
       integer, intent(out) :: ierror
     end subroutine mpi_waitsome

     llective communication

     subroutine mpi_barrier(comm, ierror) 
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
     end subroutine mpi_barrier

     subroutine mpi_op_create( function, commute, op, ierror) 
       external function 
       logical, intent(in) :: commute
       integer, intent(out) :: op
       integer, intent(out) :: ierror
     end subroutine mpi_op_create

     subroutine mpi_op_free( op, ierror) 
       integer, intent(in) :: op
       integer, intent(out) :: ierror
     end subroutine mpi_op_free

     mmunicators

     subroutine mpi_attr_delete(comm, keyval, ierror)
       integer, intent(in) :: comm
       integer, intent(in) :: keyval
       integer, intent(out) :: ierror
     end subroutine mpi_attr_delete

     subroutine mpi_attr_get(comm, keyval, attribute_val, flag, ierror) 
       integer, intent(in) :: comm
       integer, intent(in) :: keyval
       integer, intent(out) :: attribute_val
       integer, intent(out) :: ierror
       logical, intent(out) :: flag
     end subroutine mpi_attr_get

     subroutine mpi_attr_put(comm, keyval, attribute_val, ierror)
       integer, intent(in) :: comm
       integer, intent(in) :: keyval
       integer, intent(in) :: attribute_val
       integer, intent(out) :: ierror
     end subroutine mpi_attr_put

     subroutine mpi_comm_compare(comm1, comm2, result, ierror)
       integer, intent(in) :: comm1
       integer, intent(in) :: comm2
       integer, intent(out) :: result
       integer, intent(out) :: ierror
     end subroutine mpi_comm_compare

     subroutine mpi_comm_create(comm, group, newcomm, ierror)
       integer, intent(in) :: comm
       integer, intent(in) :: group
       integer, intent(out) :: newcomm
       integer, intent(out) :: ierror
     end subroutine mpi_comm_create

     subroutine mpi_comm_dup(comm, newcomm, ierror)
       integer, intent(in) :: comm
       integer, intent(out) :: newcomm
       integer, intent(out) :: ierror
     end subroutine mpi_comm_dup

     subroutine mpi_comm_free(comm, ierror)
       integer, intent(inout) :: comm
       integer, intent(out) :: ierror
     end subroutine mpi_comm_free

     subroutine mpi_comm_group(comm, group, ierror)
       integer, intent(in) :: comm
       integer, intent(out) :: group
       integer, intent(out) :: ierror
     end subroutine mpi_comm_group

     subroutine mpi_comm_rank(comm, rank, ierror)
       integer, intent(in) :: comm
       integer, intent(out) :: rank
       integer, intent(out) :: ierror
     end subroutine mpi_comm_rank

     subroutine mpi_comm_remote_group(comm, group, ierror)
       integer, intent(in) :: comm
       integer, intent(out) :: group
       integer, intent(out) :: ierror
     end subroutine mpi_comm_remote_group

     subroutine mpi_comm_remote_size(comm, size, ierror)
       integer, intent(in) :: comm
       integer, intent(out) :: size
       integer, intent(out) :: ierror
     end subroutine mpi_comm_remote_size

     subroutine mpi_comm_size(comm, size, ierror)
       integer, intent(in) :: comm
       integer, intent(out) :: size
       integer, intent(out) :: ierror
     end subroutine mpi_comm_size

     subroutine mpi_comm_split(comm, color, key, newcomm, ierror)
       integer, intent(in) :: comm
       integer, intent(in) :: color
       integer, intent(in) :: key
       integer, intent(out) :: newcomm
       integer, intent(out) :: ierror
     end subroutine mpi_comm_split

     subroutine mpi_comm_test_inter(comm, flag, ierror)
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
       logical, intent(out) :: flag
     end subroutine mpi_comm_test_inter

     subroutine mpi_group_compare(group1, group2, result, ierror)
       integer, intent(in) :: group1
       integer, intent(in) :: group2
       integer, intent(out) :: result
       integer, intent(out) :: ierror
     end subroutine mpi_group_compare

     subroutine mpi_group_difference(group1, group2, newgroup, ierror) 
       integer, intent(in) :: group1
       integer, intent(in) :: group2
       integer, intent(out) :: newgroup
       integer, intent(out) :: ierror
     end subroutine mpi_group_difference

     subroutine mpi_group_excl(group, n, ranks, newgroup, ierror)
       integer, intent(in) :: group
       integer, intent(in) :: n
       integer, intent(in) :: ranks(*)
       integer, intent(out) :: newgroup
       integer, intent(out) :: ierror
     end subroutine mpi_group_excl

     subroutine mpi_group_free(group, ierror)
       integer, intent(inout) :: group
       integer, intent(out) :: ierror
     end subroutine mpi_group_free

     subroutine mpi_group_incl(group, n, ranks, newgroup, ierror)
       integer, intent(in) :: group
       integer, intent(in) :: n
       integer, intent(in) :: ranks(*)
       integer, intent(out) :: newgroup
       integer, intent(out) :: ierror
     end subroutine mpi_group_incl

     subroutine mpi_group_intersection(group1, group2, newgroup, ierror) 
       integer, intent(in) :: group1
       integer, intent(in) :: group2
       integer, intent(out) :: newgroup
       integer, intent(out) :: ierror
     end subroutine mpi_group_intersection

     subroutine mpi_group_range_excl(group, n, ranges, newgroup, ierror) 
       integer, intent(in) :: group
       integer, intent(in) :: n
       integer, intent(in) :: ranges(3,*)
       integer, intent(out) :: newgroup
       integer, intent(out) :: ierror
     end subroutine mpi_group_range_excl

     subroutine mpi_group_range_incl(group, n, ranges, newgroup, ierror) 
       integer, intent(in) :: group
       integer, intent(in) :: n
       integer, intent(in) :: ranges(3,*)
       integer, intent(out) :: newgroup
       integer, intent(out) :: ierror
     end subroutine mpi_group_range_incl

     subroutine mpi_group_rank(group, rank, ierror)
       integer, intent(in) :: group
       integer, intent(out) :: rank
       integer, intent(out) :: ierror
     end subroutine mpi_group_rank

     subroutine mpi_group_size(group, size, ierror)
       integer, intent(in) :: group
       integer, intent(out) :: size
       integer, intent(out) :: ierror
     end subroutine mpi_group_size

     subroutine mpi_group_translate_ranks(group1, n, ranks1, &
          group2, ranks2, ierror) 
       integer, intent(in) :: group1
       integer, intent(in) :: n
       integer, intent(in) :: ranks1(*)
       integer, intent(in) :: group2
       integer, intent(out) :: ranks2(*)
       integer, intent(out) :: ierror
     end subroutine mpi_group_translate_ranks

     subroutine mpi_group_union(group1, group2, newgroup, ierror)
       integer, intent(in) :: group1
       integer, intent(in) :: group2
       integer, intent(out) :: newgroup
       integer, intent(out) :: ierror
     end subroutine mpi_group_union

     subroutine mpi_intercomm_create(local_comm, local_leader,      &
          peer_comm, remote_leader, tag, newintercomm, ierror) 
       integer, intent(in) :: local_comm
       integer, intent(in) :: local_leader
       integer, intent(in) :: peer_comm
       integer, intent(in) :: remote_leader
       integer, intent(in) :: tag
       integer, intent(out) :: newintercomm
       integer, intent(out) :: ierror
     end subroutine mpi_intercomm_create

     subroutine mpi_intercomm_merge(intercomm, high, intracomm, ierror) 
       integer, intent(in) :: intercomm
       integer, intent(out) :: intracomm
       integer, intent(out) :: ierror
       logical, intent(in) :: high
     end subroutine mpi_intercomm_merge

     subroutine mpi_keyval_create(copy_fn, delete_fn, keyval, &
          extra_state, ierror) 
       external copy_fn, delete_fn 
       integer, intent(out) :: keyval
       integer, intent(in) :: extra_state
       integer, intent(out) :: ierror
     end subroutine mpi_keyval_create

     subroutine mpi_keyval_free(keyval, ierror)
       integer, intent(inout) :: keyval
       integer, intent(out) :: ierror
     end subroutine mpi_keyval_free

     pology

     subroutine mpi_cartdim_get(comm, ndims, ierror)
       integer, intent(in) :: comm
       integer, intent(out) :: ndims
       integer, intent(out) :: ierror
     end subroutine mpi_cartdim_get

     subroutine mpi_cart_coords(comm, rank, maxdims, coords, ierror) 
       integer, intent(in) :: comm
       integer, intent(in) :: rank
       integer, intent(in) :: maxdims
       integer, intent(out) :: coords(*)
       integer, intent(out) :: ierror
     end subroutine mpi_cart_coords

     subroutine mpi_cart_create(comm_old, ndims, dims, periods,     &
          reorder, comm_cart, ierror) 
       integer, intent(in) :: comm_old
       integer, intent(in) :: ndims
       integer, intent(in) :: dims(*)
       integer, intent(out) :: comm_cart
       integer, intent(out) :: ierror
       logical, intent(in) :: periods(*)
       logical, intent(in) :: reorder
     end subroutine mpi_cart_create

     subroutine mpi_cart_get(comm, maxdims, dims, periods, coords, ierror) 
       integer, intent(in) :: comm
       integer, intent(in) :: maxdims
       integer, intent(in) :: dims(*)
       integer, intent(out) :: coords(*)
       integer, intent(out) :: ierror
       logical, intent(out) :: periods(*)
     end subroutine mpi_cart_get

     subroutine mpi_cart_map(comm, ndims, dims, periods, newrank, ierror) 
       integer, intent(in) :: comm
       integer, intent(in) :: ndims
       integer, intent(in) :: dims(*)
       integer, intent(out) :: newrank
       integer, intent(out) :: ierror
       logical, intent(in) :: periods(*)
     end subroutine mpi_cart_map

     subroutine mpi_cart_rank(comm, coords, rank, ierror)
       integer, intent(in) :: comm
       integer, intent(in) :: coords(*)
       integer, intent(out) :: rank
       integer, intent(out) :: ierror
     end subroutine mpi_cart_rank

     subroutine mpi_cart_shift(comm, direction, disp, rank_source,  &
          rank_dest, ierror) 
       integer, intent(in) :: comm
       integer, intent(in) :: direction
       integer, intent(in) :: disp
       integer, intent(out) :: rank_source
       integer, intent(out) :: rank_dest
       integer, intent(out) :: ierror
     end subroutine mpi_cart_shift

     subroutine mpi_cart_sub(comm, remain_dims, newcomm, ierror)
       integer, intent(in) :: comm
       integer, intent(out) :: newcomm
       integer, intent(out) :: ierror
       logical, intent(in) :: remain_dims(*)
     end subroutine mpi_cart_sub

     subroutine mpi_dims_create(nnodes, ndims, dims, ierror)
       integer, intent(in) :: nnodes
       integer, intent(in) :: ndims
       integer, intent(in) :: dims(*)
       integer, intent(out) :: ierror
     end subroutine mpi_dims_create

     subroutine mpi_graphdims_get(comm, nnodes, nedges, ierror)
       integer, intent(in) :: comm
       integer, intent(out) :: nnodes
       integer, intent(out) :: nedges
       integer, intent(out) :: ierror
     end subroutine mpi_graphdims_get

     subroutine mpi_graph_create(comm_old, nnodes, index, edges,    &
          reorder, comm_graph, ierror) 
       integer, intent(in) :: comm_old
       integer, intent(in) :: nnodes
       integer, intent(in) :: index(*)
       integer, intent(in) :: edges(*)
       integer, intent(out) :: comm_graph
       integer, intent(out) :: ierror
       logical, intent(in) :: reorder
     end subroutine mpi_graph_create

     subroutine mpi_graph_get(comm, maxindex, maxedges, index, edges, ierror) 
       integer, intent(in) :: comm
       integer, intent(in) :: maxindex
       integer, intent(in) :: maxedges
       integer, intent(in) :: index(*)
       integer, intent(in) :: edges(*)
       integer, intent(out) :: ierror
     end subroutine mpi_graph_get

     subroutine mpi_graph_map(comm, nnodes, index, edges, newrank, ierror) 
       integer, intent(in) :: comm
       integer, intent(in) :: nnodes
       integer, intent(in) :: index(*)
       integer, intent(in) :: edges(*)
       integer, intent(out) :: newrank
       integer, intent(out) :: ierror
     end subroutine mpi_graph_map

     subroutine mpi_graph_neighbors(comm, rank, maxneighbors, neighbors, &
          ierror) 
       integer, intent(in) :: comm
       integer, intent(in) :: rank
       integer maxneighbors
       integer, intent(out) :: neighbors(*)
       integer, intent(out) :: ierror
     end subroutine mpi_graph_neighbors

     subroutine mpi_graph_neighbors_count(comm, rank, nneighbors, ierror) 
       integer, intent(in) :: comm
       integer, intent(in) :: rank
       integer, intent(out) :: nneighbors
       integer, intent(out) :: ierror
     end subroutine mpi_graph_neighbors_count

     subroutine mpi_topo_test(comm, status, ierror)
       integer, intent(in) :: comm
       integer, intent(out) :: status
       integer, intent(out) :: ierror
     end subroutine mpi_topo_test

     vironmental inquiry

     subroutine mpi_abort(comm, errorcode, ierror)
       integer, intent(in) :: comm
       integer, intent(in) :: errorcode
       integer, intent(out) :: ierror
     end subroutine mpi_abort

     subroutine mpi_errhandler_create(function, errhandler, ierror)
       external function 
       integer, intent(out) :: errhandler
       integer, intent(out) :: ierror
     end subroutine mpi_errhandler_create

     subroutine mpi_errhandler_free(errhandler, ierror)
       integer, intent(in) :: errhandler
       integer, intent(out) :: ierror
     end subroutine mpi_errhandler_free

     subroutine mpi_errhandler_get(comm, errhandler, ierror)
       integer, intent(in) :: comm
       integer, intent(out) :: errhandler
       integer, intent(out) :: ierror
     end subroutine mpi_errhandler_get

     subroutine mpi_errhandler_set(comm, errhandler, ierror)
       integer, intent(in) :: comm
       integer, intent(in) :: errhandler
       integer, intent(out) :: ierror
     end subroutine mpi_errhandler_set

     subroutine mpi_error_class(errorcode, errorclass, ierror)
       integer, intent(in) :: errorcode
       integer, intent(out) :: errorclass
       integer, intent(out) :: ierror
     end subroutine mpi_error_class

     subroutine mpi_error_string(errorcode, string, resultlen, ierror) 
       integer, intent(in) :: errorcode
       integer, intent(out) :: resultlen
       integer, intent(out) :: ierror
       character*(*) string
     end subroutine mpi_error_string

     subroutine mpi_finalize(ierror)
       integer, intent(out) :: ierror
     end subroutine mpi_finalize

     subroutine mpi_get_processor_name( name, resultlen, ierror)
       character*(*) name
       integer resultlen,ierror
     end subroutine mpi_get_processor_name

     subroutine mpi_init(ierror)
       integer, intent(out) :: ierror
     end subroutine mpi_init

     subroutine mpi_initialized(flag, ierror)
       logical, intent(out) :: flag
       integer, intent(out) :: ierror
     end subroutine mpi_initialized

     subroutine mpi_pcontrol(level)
       integer, intent(in) :: level
     end subroutine mpi_pcontrol

     subroutine mpi_address(location, address, ierror)
       <type>, intent(in) :: location(dim1)
       integer, intent(out) :: address
       integer, intent(out) :: ierror
     end subroutine mpi_address

     subroutine mpi_bsend(buf, count, datatype, dest, tag, comm, ierror) 
       <type>, intent(in) :: buf(dim1)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
     end subroutine mpi_bsend

     subroutine mpi_bsend_init(buf, count, datatype, dest, tag,   &
          comm, request, ierror) 
       <type>, intent(in) :: buf(dim1)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
     end subroutine mpi_bsend_init

     subroutine mpi_buffer_attach( buffer, size, ierror)
       <type>, intent(in) :: buffer(dim1)
       integer, intent(in) :: size
       integer, intent(out) :: ierror
     end subroutine mpi_buffer_attach

     subroutine mpi_buffer_detach( buffer_addr, size, ierror)
       <type>, intent(out) :: buffer_addr(dim1)
       integer, intent(out) :: size
       integer, intent(out) :: ierror
     end subroutine mpi_buffer_detach

     subroutine mpi_ibsend(buf, count, datatype, dest, tag, comm, &
          request, ierror) 
       <type>, intent(in) :: buf(dim1)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
     end subroutine mpi_ibsend

     subroutine mpi_irecv(buf, count, datatype, source, tag,      &
          comm, request, ierror) 
       <type>, intent(out) :: buf(dim1)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: source
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
     end subroutine mpi_irecv

     subroutine mpi_irsend(buf, count, datatype, dest, tag, comm, &
          request, ierror) 
       <type>, intent(in) :: buf(dim1)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
     end subroutine mpi_irsend

     subroutine mpi_isend(buf, count, datatype, dest, tag, comm,  &
          request, ierror) 
       <type>, intent(in) :: buf(dim1)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
     end subroutine mpi_isend

     subroutine mpi_issend(buf, count, datatype, dest, tag, comm, &
          request, ierror) 
       <type>, intent(in) :: buf(dim1)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
     end subroutine mpi_issend

     subroutine mpi_recv(buf, count, datatype, source, tag, comm, &
          status, ierror) 
       use ModMpiOrig, only: mpi_status_size
       <type>, intent(out) :: buf(dim1)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: source
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: status(mpi_status_size)
       integer, intent(out) :: ierror
     end subroutine mpi_recv

     subroutine mpi_recv_init(buf, count, datatype, source, tag,  &
          comm, request, ierror) 
       <type>, intent(out) :: buf(dim1)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: source
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
     end subroutine mpi_recv_init

     subroutine mpi_rsend(buf, count, datatype, dest, tag, comm, ierror) 
       <type>, intent(in) :: buf(dim1)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
     end subroutine mpi_rsend

     subroutine mpi_rsend_init(buf, count, datatype, dest, tag,   &
          comm, request, ierror) 
       <type>, intent(in) :: buf(dim1)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
     end subroutine mpi_rsend_init

     subroutine mpi_send(buf, count, datatype, dest, tag, comm, ierror) 
       <type>, intent(in) :: buf(dim1)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
     end subroutine mpi_send


     subroutine mpi_sendrecv_replace(buf, count, datatype, dest,  &
          sendtag, source, recvtag, comm, status, ierror) 
       use ModMpiOrig, only: mpi_status_size
       <type>, intent(inout) :: buf(dim1)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: sendtag
       integer, intent(in) :: source
       integer, intent(in) :: recvtag
       integer, intent(in) :: comm
       integer, intent(out) :: status(mpi_status_size)
       integer, intent(out) :: ierror
     end subroutine mpi_sendrecv_replace

     subroutine mpi_send_init(buf, count, datatype, dest, tag,    &
          comm, request, ierror) 
       <type>, intent(in) :: buf(dim1)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
     end subroutine mpi_send_init

     subroutine mpi_ssend(buf, count, datatype, dest, tag, comm, ierror) 
       <type>, intent(in) :: buf(dim1)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
     end subroutine mpi_ssend

     subroutine mpi_ssend_init(buf, count, datatype, dest, tag,   &
          comm, request, ierror) 
       <type>, intent(in) :: buf(dim1)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: dest
       integer, intent(in) :: tag
       integer, intent(in) :: comm
       integer, intent(out) :: request
       integer, intent(out) :: ierror
     end subroutine mpi_ssend_init

     subroutine mpi_bcast(buffer, count, datatype, root, comm, ierror) 
       <type>, intent(inout) :: buffer(dim1)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
     end subroutine mpi_bcast

     subroutine mpi_sendrecv(sendbuf, sendcount, sendtype, dest,  &
          sendtag, recvbuf, recvcount, recvtype, source, recvtag, &
          comm, status, ierror) 
       use ModMpiOrig, only: mpi_status_size
       <type>, intent(in) :: sendbuf(dim1)
       <type>, intent(out) :: recvbuf(dim1)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: dest
       integer, intent(in) :: sendtag
       integer, intent(in) :: recvcount
       integer, intent(in) :: recvtype
       integer, intent(in) :: source
       integer, intent(in) :: recvtag
       integer, intent(in) :: comm
       integer, intent(out) :: status(mpi_status_size)
       integer, intent(out) :: ierror
     end subroutine mpi_sendrecv

     subroutine mpi_allgather(sendbuf, sendcount, sendtype,       &
          recvbuf, recvcount, recvtype, comm, ierror) 
       <type>, intent(in) :: sendbuf(dim1)
       <type>, intent(out) :: recvbuf(dim2)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcount
       integer, intent(in) :: recvtype
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
     end subroutine mpi_allgather

     subroutine mpi_allgatherv(sendbuf, sendcount, sendtype,      &
          recvbuf, recvcounts, displs, recvtype, comm, ierror) 
       <type>, intent(in) :: sendbuf(dim1)
       <type>, intent(out) :: recvbuf(dim1)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcounts(:)
       integer, intent(in) :: displs(:)
       integer, intent(in) :: recvtype
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
     end subroutine mpi_allgatherv

     subroutine mpi_allreduce(sendbuf, recvbuf, count, datatype,  &
          op, comm, ierror) 
       <type>, intent(in) :: sendbuf(dim1)
       <type>, intent(out) :: recvbuf(dim1)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: op
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
     end subroutine mpi_allreduce

     subroutine mpi_alltoall(sendbuf, sendcount, sendtype,        &
          recvbuf, recvcount, recvtype, comm, ierror) 
       <type>, intent(in) :: sendbuf(dim1)
       <type>, intent(out) :: recvbuf(dim2)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcount
       integer, intent(in) :: recvtype
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
     end subroutine mpi_alltoall

     subroutine mpi_alltoallv(sendbuf, sendcounts, sdispls,       &
          sendtype, recvbuf, recvcounts, rdispls, recvtype, comm, &
          ierror) 
       <type>, intent(in) :: sendbuf(dim1)
       <type>, intent(out) :: recvbuf(dim1)
       integer, intent(in) :: sendcounts(:)
       integer, intent(in) :: sdispls(:)
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcounts(:)
       integer, intent(in) :: rdispls(:)
       integer, intent(in) :: recvtype
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
     end subroutine mpi_alltoallv

     subroutine mpi_gather(sendbuf, sendcount, sendtype, recvbuf, &
          recvcount, recvtype, root, comm, ierror) 
       <type>, intent(in) :: sendbuf(dim1)
       <type>, intent(out) :: recvbuf(dim2)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcount
       integer, intent(in) :: recvtype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
     end subroutine mpi_gather

     subroutine mpi_gatherv(sendbuf, sendcount, sendtype,         &
          recvbuf, recvcounts, displs, recvtype, root, comm, ierror) 
       <type>, intent(in) :: sendbuf(dim1)
       <type>, intent(out) :: recvbuf(dim1)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcounts(:)
       integer, intent(in) :: displs(:)
       integer, intent(in) :: recvtype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
     end subroutine mpi_gatherv

     subroutine mpi_reduce(sendbuf, recvbuf, count, datatype, op, &
          root, comm, ierror) 
       <type>, intent(in) :: sendbuf(dim1)
       <type>, intent(out) :: recvbuf(dim1)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: op
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
     end subroutine mpi_reduce

     subroutine mpi_reduce_scatter(sendbuf, recvbuf, recvcounts,  &
          datatype, op, comm, ierror) 
       <type>, intent(in) :: sendbuf(dim1)
       <type>, intent(out) :: recvbuf(dim1)
       integer, intent(in) :: recvcounts(:)
       integer, intent(in) :: datatype
       integer, intent(in) :: op
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
     end subroutine mpi_reduce_scatter

     subroutine mpi_scan(sendbuf, recvbuf, count, datatype, op,   &
          comm, ierror) 
       <type>, intent(in) :: sendbuf(dim1)
       <type>, intent(out) :: recvbuf(dim2)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: op
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
     end subroutine mpi_scan

     subroutine mpi_scatter(sendbuf, sendcount, sendtype,         &
          recvbuf, recvcount, recvtype, root, comm, ierror) 
       <type>, intent(in) :: sendbuf(dim1)
       <type>, intent(out) :: recvbuf(dim2)
       integer, intent(in) :: sendcount
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcount
       integer, intent(in) :: recvtype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
     end subroutine mpi_scatter

     subroutine mpi_scatterv(sendbuf, sendcounts, displs,         &
          sendtype, recvbuf, recvcount, recvtype, root, comm, ierror) 
       <type>, intent(in) :: sendbuf(dim1)
       <type>, intent(out) :: recvbuf(dim2)
       integer, intent(in) :: sendcounts(:)
       integer, intent(in) :: displs(:)
       integer, intent(in) :: sendtype
       integer, intent(in) :: recvcount
       integer, intent(in) :: recvtype
       integer, intent(in) :: root
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
     end subroutine mpi_scatterv

     subroutine mpi_pack(inbuf, incount, datatype, outbuf, outsize, &
          position, comm, ierror)
       <type>, intent(in) :: inbuf(dim1)
       <type>, intent(out) :: outbuf(dim2)
       integer, intent(in) :: incount
       integer, intent(in) :: datatype
       integer, intent(in) :: outsize
       integer, intent(inout) :: position
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
     end subroutine mpi_pack

     subroutine mpi_unpack(inbuf, insize, position, outbuf,       &
          outcount, datatype, comm, ierror)
       <type>, intent(in) :: inbuf(dim1)
       <type>, intent(out) :: outbuf(dim2)
       integer, intent(in) :: insize
       integer, intent(inout) :: position
       integer, intent(in) :: outcount
       integer, intent(in) :: datatype
       integer, intent(in) :: comm
       integer, intent(out) :: ierror
     end subroutine mpi_unpack

     subroutine mpi_file_read(fh, buf, count, datatype, status, ierror)
       use ModMpiOrig, only: mpi_status_size
       integer, intent(in) :: fh
       <type>,  intent(out):: buf(dim1)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: status(mpi_status_size)
       integer, intent(out) :: ierror
     end subroutine mpi_file_read

     subroutine mpi_file_write(fh, buf, count, datatype, status, ierror)
       use ModMpiOrig, only: mpi_status_size
       integer, intent(in) :: fh
       <type>,  intent(in) :: buf(dim1)
       integer, intent(in) :: count
       integer, intent(in) :: datatype
       integer, intent(in) :: status(mpi_status_size)
       integer, intent(out) :: ierror
     end subroutine mpi_file_write

  end interface            ! These two lines are here so that 
end module ModMpiTemplate  ! EMACS can indent the code properly
