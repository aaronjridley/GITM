!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!-*- F90 -*- so emacs thinks this is an f90 file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! NAME: AB_COMM_module_mpi
!
! PURPOSE: a module which implements the communication neccessary for 
!          AB2D using mpi.
!
! HISTORY:
!  10/20/99 Robert Oehmke: created
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module AB_COMM_module
  use AB_ARRAY_module
  use AB_ERROR_module
  implicit none

  type :: AB_COMM
     integer :: context
     integer :: me,np,max_pn
  end type AB_COMM

  type :: AB_COMM_XCHNG
     type(AB_COMM),pointer :: comm
     integer, pointer :: snd_req(:)
     integer, pointer :: rcv_req(:)
     integer, pointer :: status(:,:)
     integer :: tag
  end type AB_COMM_XCHNG

contains


subroutine AB_COMM_create(c,context)
  use ModMpi
  type(AB_COMM), intent(out) ::c
  integer, intent(in) :: context
  integer ierror


  c%context=context
  call MPI_COMM_RANK(context, c%me, ierror)
  call MPI_COMM_SIZE(context, c%np, ierror)
  c%max_pn=c%np-1

end subroutine AB_COMM_create

subroutine AB_COMM_get_num_procs(c,num_procs)
  type(AB_COMM), intent(in) ::c
  integer, intent(out) :: num_procs

  num_procs=c%np

end subroutine AB_COMM_get_num_procs


subroutine AB_COMM_get_my_proc_num(c,me)
  type(AB_COMM), intent(in) ::c
  integer, intent(out) :: me

  me=c%me

end subroutine AB_COMM_get_my_proc_num


subroutine AB_COMM_XCHNG_create(xchng,c,tag,ok)
  use ModMpi
  type(AB_COMM), intent(in),target ::c
  integer, intent(in) :: tag
  type(AB_COMM_XCHNG), intent(out) :: xchng
  logical, intent(out) :: ok
  integer i,max_pn, ierror

  ! intialize status flag
  ok=.true.

  ! to make things easier get some variables
  max_pn=c%max_pn

  ! initialize variables
  xchng%comm=>c
  xchng%tag=tag

  ! allocate global status array used in transfer procedures
  allocate(xchng%status(MPI_STATUS_SIZE,0:max_pn), stat=ierror)
  if (ierror .ne. 0) then
     ok=.false.
     call AB_ERROR_set("AB_COMM_XCHNG_create","allocate error ",ierror)
     return
  endif

  ! allocate global receive request array used in transfer procedures
  allocate(xchng%rcv_req(0:max_pn), stat=ierror)
  if (ierror .ne. 0) then
     ok=.false.
     call AB_ERROR_set("AB_COMM_XCHNG_create","allocate error ",ierror)
     return
  endif
  xchng%rcv_req=MPI_REQUEST_NULL

  ! allocate global send request array used in transfer procedures
  allocate(xchng%snd_req(0:max_pn), stat=ierror)
  if (ierror .ne. 0) then
     ok=.false.
     call AB_ERROR_set("AB_COMM_XCHNG_create","allocate error ",ierror)
     return
  endif
  xchng%snd_req=MPI_REQUEST_NULL

end subroutine AB_COMM_XCHNG_create

subroutine AB_COMM_XCHNG_REAL_start(xchng,snd_size,snd, &
                                    rcv_size,rcv, &
                                    ok)
  use ModMpi
  type(AB_COMM_XCHNG), intent(inout) :: xchng
  integer, dimension(0:) :: snd_size, rcv_size
  type(AB_ARRAY_REAL), dimension(0:) :: snd, rcv 
  logical, intent(out) :: ok
  integer :: p, buf_size, ierror
  integer :: max_pn,me,context,tag,np

  ! intialize status flag
  ok=.true.

  ! get some variables to make things clearer
  max_pn=xchng%comm%max_pn
  np=xchng%comm%np
  me=xchng%comm%me
  context=xchng%comm%context
  tag=xchng%tag

  ! make sure there are no outstanding receive requests
  call MPI_waitall(np,xchng%rcv_req,xchng%status,ierror)
  if (ierror .ne. MPI_SUCCESS) then
     ok=.false.
     call AB_ERROR_set("AB_COMM_XCHNG_REAL_start","MPI error ",ierror)
     return
  endif

  ! initialize receive requests to null
  xchng%rcv_req=MPI_REQUEST_NULL

  ! setup buffers to receive
  do p=0,max_pn
     if (p .ne. me) then
        buf_size=rcv_size(p)
        if (buf_size .ne. 0) then
           call MPI_irecv(rcv(p)%array(1), buf_size, MPI_REAL, &
                p, tag, context, &
                xchng%rcv_req(p), ierror) 
           if (ierror .ne. MPI_SUCCESS) then
              ok=.false.
              call AB_ERROR_set("AB_COMM_XCHNG_REAL_start","MPI error ",ierror)
              return
           endif
        endif
     endif
  enddo

  ! make sure there are no outstanding send requests
  call MPI_waitall(np,xchng%snd_req,xchng%status,ierror)
  if (ierror .ne. MPI_SUCCESS) then
     ok=.false.
     call AB_ERROR_set("AB_COMM_XCHNG_REAL_start","MPI error ",ierror)
     return
  endif

  ! initialize send requests to null
  xchng%snd_req=MPI_REQUEST_NULL
  
  ! send data
  do p=0,max_pn
     if (p .ne. me) then
        buf_size=snd_size(p)
        if (buf_size .ne. 0) then
           call MPI_isend(snd(p)%array(1), buf_size, MPI_REAL, &
                          p, tag, context, &
                          xchng%snd_req(p), ierror)                     
           if (ierror .ne. MPI_SUCCESS) then
              ok=.false.
              call AB_ERROR_set("AB_COMM_XCHNG_REAL_start","MPI error ",ierror)
              return
           endif
        endif
     endif
  enddo


end subroutine AB_COMM_XCHNG_REAL_start

subroutine AB_COMM_XCHNG_INT_2D_start(xchng,snd_size,snd_size_mult,snd, &
                                    rcv_size,rcv_size_mult,rcv, &
                                    ok)
  use ModMpi
  type(AB_COMM_XCHNG), intent(inout) :: xchng
  integer, dimension(0:) :: snd_size, rcv_size
  type(AB_ARRAY_INT_2D), dimension(0:) :: snd, rcv
  integer, intent(in) :: rcv_size_mult, snd_size_mult 
  logical, intent(out) :: ok
  integer :: p, buf_size, ierror
  integer :: max_pn,me,context,tag,np

  ! intialize status flag
  ok=.true.

  ! get some variables to make things clearer
  max_pn=xchng%comm%max_pn
  np=xchng%comm%np
  me=xchng%comm%me
  context=xchng%comm%context
  tag=xchng%tag

  ! make sure there are no outstanding receive requests
  call MPI_waitall(np,xchng%rcv_req,xchng%status,ierror)
  if (ierror .ne. MPI_SUCCESS) then
     ok=.false.
     call AB_ERROR_set("AB_COMM_XCHNG_INT_2D_start","MPI error ",ierror)
     return
  endif

  ! initialize receive requests to null
  xchng%rcv_req=MPI_REQUEST_NULL

  ! setup buffers to receive
  do p=0,max_pn
     if (p .ne. me) then
        buf_size=rcv_size(p)*rcv_size_mult
        if (buf_size .ne. 0) then
           call MPI_irecv(rcv(p)%array(1,1), buf_size, MPI_INTEGER, &
                p, tag, context, &
                xchng%rcv_req(p), ierror) 
           if (ierror .ne. MPI_SUCCESS) then
              ok=.false.
            call AB_ERROR_set("AB_COMM_XCHNG_INT_2D_start","MPI error ",ierror)
              return
           endif
        endif
     endif
  enddo

  ! make sure there are no outstanding send requests
  call MPI_waitall(np,xchng%snd_req,xchng%status,ierror)
  if (ierror .ne. MPI_SUCCESS) then
     ok=.false.
     call AB_ERROR_set("AB_COMM_XCHNG_INT_2D_start","MPI error ",ierror)
     return
  endif

  ! initialize send requests to null
  xchng%snd_req=MPI_REQUEST_NULL
  
  ! send data
  do p=0,max_pn
     if (p .ne. me) then
        buf_size=snd_size(p)*snd_size_mult
        if (buf_size .ne. 0) then
           call MPI_isend(snd(p)%array(1,1), buf_size, MPI_INTEGER, &
                          p, tag, context, &
                          xchng%snd_req(p), ierror)                     
           if (ierror .ne. MPI_SUCCESS) then
              ok=.false.
              call AB_ERROR_set("AB_COMM_XCHNG_INT_2D_start","MPI error ",ierror)
              return
           endif
        endif
     endif
  enddo


end subroutine AB_COMM_XCHNG_INT_2D_start

subroutine AB_COMM_XCHNG_finish_rcv(xchng,ok)
  use ModMpi
  type(AB_COMM_XCHNG), intent(inout) :: xchng
  logical, intent(out) :: ok
  integer :: ierror

  ! intialize status flag
  ok=.true.

  ! wait for all the receives to finish
  call MPI_waitall(xchng%comm%np,xchng%rcv_req,xchng%status,ierror)
  if (ierror .ne. MPI_SUCCESS) then
     ok=.false.
     call AB_ERROR_set("AB_COMM_XCHNG_finish_rcv","MPI error ",ierror)
     return
  endif

end subroutine AB_COMM_XCHNG_finish_rcv

subroutine AB_COMM_XCHNG_finish_snd(xchng,ok)
  use ModMpi
  type(AB_COMM_XCHNG), intent(inout) :: xchng
  logical, intent(out) :: ok
  integer :: ierror

  ! intialize status flag
  ok=.true.

  ! wait for all the receives to finish
  call MPI_waitall(xchng%comm%np,xchng%snd_req,xchng%status,ierror)
  if (ierror .ne. MPI_SUCCESS) then
     ok=.false.
     call AB_ERROR_set("AB_COMM_XCHNG_finish_snd","MPI error ",ierror)
     return
  endif

end subroutine AB_COMM_XCHNG_finish_snd

subroutine AB_COMM_XCHNG_destroy(xchng,ok)
  use ModMpi
  type(AB_COMM_XCHNG), intent(inout) :: xchng
  logical, intent(out) :: ok
  integer i,np,ierror

  ! intialize status flag
  ok=.true.

  ! to make things easier get some variables
  np=xchng%comm%np

  ! wait for all the sends to complete
  call MPI_waitall(np,xchng%snd_req,xchng%status,ierror)
  if (ierror .ne. MPI_SUCCESS) then
     ok=.false.
     call AB_ERROR_set("AB_COMM_XCHNG_destroy","MPI error ",ierror)
     return
  endif

  ! deallocate send request array
  deallocate(xchng%snd_req)

  ! wait for all the rcvs to complete
  call MPI_waitall(np,xchng%rcv_req,xchng%status,ierror)
  if (ierror .ne. MPI_SUCCESS) then
     ok=.false.
     call AB_ERROR_set("AB_COMM_XCHNG_destroy","MPI error ",ierror)
     return
  endif

  ! deallocate receive request array
  deallocate(xchng%rcv_req)

  ! deallocate status array
  deallocate(xchng%status)

end subroutine AB_COMM_XCHNG_destroy


end module AB_COMM_module


