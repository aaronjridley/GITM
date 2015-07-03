!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!-*- F90 -*- so emacs thinks this is an f90 file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! NAME: AB_MSG_module
!
! PURPOSE: an internal module which implements a message exchange system
!
! HISTORY:
!  2/14/01 Robert Oehmke: created
!
! NOTE:
!  For internal use only, not sturdily enough built for general use
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module AB_MSG_module
  use AB_ARRAY_module
  use AB_ERROR_module
  use AB_COMM_module
  implicit none

  integer, private, parameter :: msg_len = 4

  type private :: AB_MSG
     integer, pointer :: curr_out(:)
     type(AB_ARRAY_INT) ::out(:)
     type(AB_COMM_XCHNG) :: xchng
     integer :: me,np,max_pn
  end type AB_MSG

   interface AB_MSG_add
      module procedure AB_MSG_add7
   end interface


contains


subroutine AB_MSG_create(msg,comm,max_len,ok)
  type(AB_MSG), intent(out) :: msg
  type(AB_COMM),intent(in) :: comm
  integer, dimension(0:) :: max_len
  logical, intent(out) :: ok
  logical :: tmp_ok
  integer :: ierror
  integer :: p,max_pn


  ! initialize error return
  ok=.true.

  ! get parallel info
  call AB_COMM_get_my_proc_num(comm, msg%me)
  call AB_COMM_get_num_procs(comm, msg%np)
  msg%max_pn=msg%np-1  
  max_pn=msg%max_pn

  ! create xchng structure
  call AB_COMM_XCHNG_create(msg%xchng,comm,99,tmp_ok)
  if (.not. tmp_ok) then
     ok=.false.
     return
  endif

  ! allocate and initialize current message position
  allocate(msg%curr_out(0:max_pn), stat=ierror)
  if (ierror .ne. 0) then
     ok=.false.
     call AB_ERROR_set("AB_MSG_create","allocate error ",ierror)
     return
  endif
  msg%curr_out=0

  ! allocate outgoing message buffers
  allocate(msg%out(0:max_pn), stat=ierror)
  if (ierror .ne. 0) then
     ok=.false.
     call AB_ERROR_set("AB_MSG_create","allocate error ",ierror)
     return
  endif

  do p=0,max_pn
     allocate(msg%out(p)%array(max_len(p)), stat=ierror)
     if (ierror .ne. 0) then
        ok=.false.
        call AB_ERROR_set("AB_MSG_create","allocate error ",ierror)
        return
     endif
  enddo

end subroutine AB_MSG_create


! should add error checking to this so we don't get more
! msgs than we have space for
subroutine AB_MSG_add7(msg,cmd,to_p,v1,v2,v3,v4,v5,v6,v7)
  type(AB_MSG), intent(inout) :: msg
  integer, intent(in) :: to_p,cmd,v1,v2,v3,v4,v5,v6,v7
  integer :: pos

  pos=msg%cur_out(to_p)
  msg%out(to_p)%array(pos)=cmd
  msg%out(to_p)%array(pos+1)=v1
  msg%out(to_p)%array(pos+2)=v2
  msg%out(to_p)%array(pos+3)=v3
  msg%out(to_p)%array(pos+4)=v4
  msg%out(to_p)%array(pos+5)=v5
  msg%out(to_p)%array(pos+6)=v6
  msg%out(to_p)%array(pos+7)=v7
  msg%cur_out(to_p)=pos+8

end subroutine AB_MSG_add7



subroutine AB_MSG_destroy(msg,ok)
  type(AB_MSG), intent(inout) :: msg
  logical, intent(out) :: ok
  logical :: tmp_ok
  integer :: p


  ! initialize status
  ok=.true.

  ! destroy the exchange structure
  call AB_COMM_XCHNG_destroy(ab_xfer%xchng,tmp_ok)
  if (.not. tmp_ok) then
     ok=.false.
     return
  endif

 ! deallocate buffers
 do p=0,msg%max_pn
    deallocate(msg%out(p)%array)
 enddo

 ! deallocate other memory
 deallocate(msg%curr_out)

end subroutine AB_MSG_destroy





end module AB_COMM_module


