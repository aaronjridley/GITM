!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!-*- F90 -*- so emacs thinks this is an f90 file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! NAME: AB_XFER_module
!
! PURPOSE: a module which implements the user communication for 
!          AB.
!
! HISTORY:
!  1/05/01 Robert Oehmke: created
!
! NOTES: 
!  Someday add function which directly copies from one block to another
!  on same proc. to reduce the amount of buffer space we need. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module AB_XFER_module
  use AB_module
  use AB_COMM_module
  use AB_ARRAY_module
  use AB_ERROR_module
  implicit none

  type AB_XFER
     private
     integer :: size(ab_num_nbrs)
     integer :: tag
     integer :: me,np,max_pn
     integer, pointer :: pos(:)
     type(AB_COMM), pointer :: comm 
     type(AB_COMM_XCHNG) :: xchng
     type(AB_ARRAY_REAL), pointer :: snd_buf(:), rcv_buf(:) 
     integer, pointer :: snd_buf_size(:), rcv_buf_size(:)
     type (AB_GRP),pointer :: grp
     integer, pointer :: num_nbrs(:) ! indexed beg. at 0
     type(AB_ARRAY_INT_2D), pointer :: nbr_order(:) ! indexed beg. at 0 
                                                     ! (1,x) -> blk number to
                                                     ! (2,x) -> dir to
                                                     ! (3,x) -> dir from
  end type AB_XFER

contains


  ! setup msg aggeragation info needed by xfers
  subroutine setup_for_xfers(xfer, ok)
    type (AB_XFER), intent(inout) :: xfer
    logical, intent(inout) :: ok
    integer :: ierror,i,j,tmp_proc, pos, buf_size
    type(AB_ARRAY_INT_2D), pointer :: snd_bufs(:) 
    type (AB_NBR_ITER) :: iter
    logical :: tmp_ok
    integer :: ind,dir,nbr_p,nbr_i, nbr_dir
    logical :: done, pole

    ! initialize stat
    ok=.true.

    ! create the neighbor iterator
    call AB_NBR_ITER_create(iter,xfer%grp)

    ! Allocate space for the number of neighboring blocks belonging to each 
    ! processor
    allocate(xfer%num_nbrs(0:xfer%max_pn), stat=ierror)
    if (ierror .ne. 0) then
       ok=.false.
       call AB_ERROR_set("setup_for_xfers","allocate error ",ierror)
       return
    endif

    ! Calculate the number of neighbors belonging to each processor
    xfer%num_nbrs=0
    call AB_NBR_ITER_reset(iter,ind,dir,nbr_p,nbr_i,done)
    do while (.not. done)
       if (nbr_p .ne. -1) then
          xfer%num_nbrs(nbr_p)=xfer%num_nbrs(nbr_p)+1
       endif

       call AB_NBR_ITER_next(iter,ind,dir,nbr_p,nbr_i,done)
    enddo

    ! Allocate space for nbr order info
    allocate(xfer%nbr_order(0:xfer%max_pn), stat=ierror)
    if (ierror .ne. 0) then
       ok=.false.
       call AB_ERROR_set("setup_for_xfers","allocate error ",ierror)
       return
    endif

    ! Allocate space for nbr order snd buffers
    allocate(snd_bufs(0:xfer%max_pn), stat=ierror)
    if (ierror .ne. 0) then
       ok=.false.
       call AB_ERROR_set("setup_for_xfers","allocate error ",ierror)
       return
    endif

    ! allocate order info space and snd buffers 
    do i=0,xfer%max_pn
       allocate(xfer%nbr_order(i)%array(3,xfer%num_nbrs(i)), &
                stat=ierror)
       if (ierror .ne. 0) then
          ok=.false.
          call AB_ERROR_set("setup_for_xfers","allocate error ",ierror)
          return
       endif

       allocate(snd_bufs(i)%array(3,xfer%num_nbrs(i)), &
                stat=ierror)
       if (ierror .ne. 0) then
          ok=.false.
          call AB_ERROR_set("setup_for_xfers","allocate error ",ierror)
          return
       endif
    enddo

    ! fill snd bufs with order info
    xfer%pos=1 ! intialize all entries in xfer%pos to 1
    call AB_NBR_ITER_reset(iter,ind,dir,nbr_p,nbr_i,pole,nbr_dir,done)
    do while (.not. done)
       if (nbr_p .ne. -1) then
          pos=xfer%pos(nbr_p)
          snd_bufs(nbr_p)%array(1,pos)=nbr_i
          snd_bufs(nbr_p)%array(2,pos)=nbr_dir
          snd_bufs(nbr_p)%array(3,pos)=dir
          pos=pos+1
          xfer%pos(nbr_p)=pos
       endif

       call AB_NBR_ITER_next(iter,ind,dir,nbr_p,nbr_i,pole,nbr_dir,done)
    enddo

    ! start exchange of data
    call AB_COMM_XCHNG_INT_2D_start(xfer%xchng,xfer%num_nbrs,3,snd_bufs, &
                                    xfer%num_nbrs,3,xfer%nbr_order, &
                                    tmp_ok)
    if (.not. tmp_ok) then
       ok=.false.
       return
    endif

    ! copy same proc info
    if (xfer%num_nbrs(xfer%me) .ne. 0) then
       xfer%nbr_order(xfer%me)%array=snd_bufs(xfer%me)%array
    endif



    ! wait until we've sent everything before deallocating the send bufs
    call AB_COMM_XCHNG_finish_snd(xfer%xchng,tmp_ok)
    if (.not. tmp_ok) then
       ok=.false.
       return
    endif


    ! Deallocate send buffers
    do i=0,xfer%max_pn
       deallocate(snd_bufs(i)%array)
    enddo
    deallocate(snd_bufs)

    ! wait until we've received everything before exiting
    call AB_COMM_XCHNG_finish_rcv(xfer%xchng,tmp_ok)
    if (.not. tmp_ok) then
       ok=.false.
       return
    endif

  end subroutine setup_for_xfers


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! NAME: AB_XFER_create
 !
 ! PURPOSE: create a transfer structure. Allocate buffers, initialize 
 !          data structs, etc. 
 !
 ! INPUTS: grp    - the adaptive block group on which the transfer is to
 !                     be performed.
 !         size(ab_num_nbrs)  
 !                    - An array of sizes indexed by d where size(d) is the
 !                      number of reals being sent to the neighbor in direction
 !                      d.
 !         tag       - tag unique amongst all other concurrent MPI 
 !                     communication. Used to identify this transfer. 
 !
 ! OUTPUTS: xfer   - the newly created transfer structure.
 !          ok     - status 
 !
 ! HISTORY:
 !  10/20/99 Robert Oehmke: created
 !  12/23/99 Robert Oehmke: changed so all communication between
 !                          2 procs is aggeragated into one big message.
 !   5/05/00 Robert Oehmke: allowed msgs in different directions to vary
 !                          in size.
 !   5/11/00 Robert Oehmke: added corner communication.
 !   2/06/01 Robert Oehmke: moved to seperate module
 !
 ! NOTES:
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine AB_XFER_create(xfer, grp, size, tag, ok)
    type (AB_GRP), intent(in),target :: grp
    type (AB_XFER), intent(out) :: xfer
    integer, intent(in) :: size(ab_num_nbrs), tag
    logical, intent(out) :: ok
    logical :: tmp_ok
    integer :: ierror,i,j
    integer :: buf_size,p,num_blks,num_rcv, rcv_buf_size
    type (AB_NBR_ITER) :: iter
    integer, pointer :: rcv_order(:,:)
    integer :: ind,dir,nbr_p,nbr_i
    logical :: done

    ! intialize status 
    ok=.true.

    ! initialize variables
    xfer%grp=>grp
    xfer%size=size
    xfer%tag=tag
    xfer%comm=>g_comm ! do something better than this when you have time
    call AB_COMM_get_my_proc_num(xfer%comm, xfer%me)
    call AB_COMM_get_num_procs(xfer%comm, xfer%np)
    xfer%max_pn=xfer%np-1

    ! create the neighbor iterator
    call AB_NBR_ITER_create(iter,grp)

    ! create xchng structure
    call AB_COMM_XCHNG_create(xfer%xchng,xfer%comm,tag,tmp_ok)
    if (.not. tmp_ok) then
       ok=.false.
       return
    endif

    ! allocate pos array used in transfer procedures
    allocate(xfer%pos(0:xfer%max_pn), stat=ierror)
    if (ierror .ne. 0) then
       ok=.false.
       call AB_ERROR_set("AB_XFER_create","allocate error ",ierror)
       return
    endif

    ! Setup for transfers
    call setup_for_xfers(xfer, tmp_ok)
    if (.not. tmp_ok) then
       ok=.false.
       return
    endif

    ! setup snd buffers
    !! allocate snd buffer storage
    allocate(xfer%snd_buf(0:xfer%max_pn), stat=ierror)
    if (ierror .ne. 0) then
       ok=.false.
       call AB_ERROR_set("AB_XFER_create","allocate error ",ierror)
       return
    endif

    !! allocate array to hold size of snd buffers
    allocate(xfer%snd_buf_size(0:xfer%max_pn), stat=ierror)
    if (ierror .ne. 0) then
       ok=.false.
       call AB_ERROR_set("AB_XFER_create","allocate error ",ierror)
       return
    endif

    !! calculate size of snd buffers
    xfer%snd_buf_size=0
    call AB_NBR_ITER_reset(iter,ind,dir,nbr_p,nbr_i,done)
    do while (.not. done)
       if (nbr_p .ne. -1) then
          xfer%snd_buf_size(nbr_p)=xfer%snd_buf_size(nbr_p)+xfer%size(dir)
       endif

       call AB_NBR_ITER_next(iter,ind,dir,nbr_p,nbr_i,done)
    enddo

    !! allocate snd buffers
    do p=0,xfer%max_pn
       !!! get the size of p's buffer
       buf_size=xfer%snd_buf_size(p)

       !!! if there's a buffer to allocate than do so
       if (buf_size .ne. 0) then
          allocate(xfer%snd_buf(p)%array(buf_size),stat=ierror)
          if (ierror .ne. 0) then
             ok=.false.
             call AB_ERROR_set("AB_XFER_create","allocate error ",ierror)
             return
          endif
       else
          nullify(xfer%snd_buf(p)%array)         
       endif
    enddo

    ! setup rcv buffers
    !! allocate rcv buffer storage
    allocate(xfer%rcv_buf(0:xfer%max_pn), stat=ierror)
    if (ierror .ne. 0) then
       ok=.false.
       call AB_ERROR_set("AB_XFER_create","allocate error ",ierror)
       return
    endif

    !! allocate array to hold size of buffers
    allocate(xfer%rcv_buf_size(0:xfer%max_pn), stat=ierror)
    if (ierror .ne. 0) then
       ok=.false.
       call AB_ERROR_set("AB_XFER_create","allocate error ",ierror)
       return
    endif

    !! loop through blocks totalling the size coming into each proc
    xfer%rcv_buf_size=0
    do p=0,xfer%max_pn
       num_rcv=xfer%num_nbrs(p)
       rcv_order=>xfer%nbr_order(p)%array
       rcv_buf_size=0
       do i=1,num_rcv
          rcv_buf_size=rcv_buf_size+xfer%size(rcv_order(3,i)) 
       enddo
       xfer%rcv_buf_size(p)=rcv_buf_size 
    enddo

    !! allocate buffers
    do p=0,xfer%max_pn
       !!! get the size of p's buffer
       buf_size=xfer%rcv_buf_size(p)

       !!! if there's a buffer to allocate than do so
       if (buf_size .ne. 0) then
          allocate(xfer%rcv_buf(p)%array(buf_size),stat=ierror)
          if (ierror .ne. 0) then
             ok=.false.
             call AB_ERROR_set("AB_XFER_create","allocate error ",ierror)
             return
          endif
       else
          nullify(xfer%rcv_buf(p)%array)         
       endif
    enddo

  end subroutine AB_XFER_create

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! NAME: AB_XFER_start
 !
 ! PURPOSE: start a real transfer. 
 !
 ! INPUTS: xfer - transfer structure for which to do transfer.
 !         get     - procedure of the form get(ind,dir,pole,oarray)
 !                   get's function is to fill oarray with the
 !                   data from the block at index ind to be sent
 !                   in the neighbor in direction dir. pole is true
 !                   if the data will cross the pole.  
 !         get4    - procedure of the form 
 !                   get4(ind,dir,oarray1,oarray2,oarray3,oarray4)
 !                   CURRENTLY UNUSED.
 !
 ! OUTPUTS: ok - status
 !      
 !
 ! HISTORY:
 !  10/20/99 Robert Oehmke: created
 !  12/23/99 Robert Oehmke: changed so all communication between
 !                          2 procs is aggeragated into one big message.
 !   5/05/00 Robert Oehmke: allowed msgs in different directions to vary
 !                          in size.
 !   5/11/00 Robert Oehmke: added corner communication.
 !   2/06/01 Robert Oehmke: moved to seperate module
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine AB_XFER_start(xfer,get,get4,ok)
    type (AB_XFER), intent(inout) :: xfer
    external get4
    interface
       subroutine get(index,dir,pole,out_array)
         integer, intent(in) :: index,dir
         logical, intent(in) :: pole
         real, dimension(:),intent(out) :: out_array
       end subroutine get
    end interface
    logical, intent(out) :: ok
    logical :: tmp_ok
    integer :: ierror
    integer :: buf_size, min_buf, max_buf
    integer :: ind,dir,nbr_p,nbr_i, nbr_dir
    logical :: pole, done
    type (AB_NBR_ITER) :: iter


    ! intialize status 
    ok=.true.

    ! create a neighbor iterator
    call AB_NBR_ITER_create(iter,xfer%grp)

    ! get data from user
    xfer%pos=1 ! intialize all entries in xfer%pos to 1
    call AB_NBR_ITER_reset(iter,ind,dir,nbr_p,nbr_i,pole,nbr_dir,done)
    do while (.not. done)

       if (nbr_p .ne. -1) then
          min_buf=xfer%pos(nbr_p)
          max_buf=min_buf+xfer%size(dir)-1
          call get(ind,dir,pole,xfer%snd_buf(nbr_p)%array(min_buf:max_buf))
          xfer%pos(nbr_p)=max_buf+1
       endif

       call AB_NBR_ITER_next(iter,ind,dir,nbr_p,nbr_i,pole,nbr_dir,done)
    enddo

    ! start the exchange of data
    call AB_COMM_XCHNG_REAL_start(xfer%xchng, &
                                  xfer%snd_buf_size,xfer%snd_buf, &
                                  xfer%rcv_buf_size,xfer%rcv_buf, &
                                  tmp_ok)
    if (.not. tmp_ok) then
       ok=.false.
       return
    endif

  end subroutine AB_XFER_start


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! NAME: AB_XFER_finish
 !
 ! PURPOSE: wait for a real transfer to finish and put the transfered data 
 !          back where the user wants it.           
 !
 ! INPUTS: xfer - transfer structure for which to finish transfer.
 !         put     - procedure of the form put(ind,dir,iarray)
 !                   put's function is to take the data in iarray 
 !                   coming from the neighbor in direction dir and
 !                   put it into the block at index ind.
 !                   in the neighbor in direction dir.   
 !         put4    - procedure of the form 
 !                   put4(ind,dir,iarray1,iarray2,iarray3,iarray4)
 !                   currently unused.
 ! 
 ! OUTPUTS: ok -status
 !
 ! HISTORY:
 !  10/20/99 Robert Oehmke: created
 !  12/23/99 Robert Oehmke: changed so all communication between
 !                          2 procs is aggeragated into one big message.
 !   5/05/00 Robert Oehmke: allowed msgs in different directions to vary
 !                          in size.
 !   5/11/00 Robert Oehmke: added corner communication.
 !   2/06/01 Robert Oehmke: moved to seperate module
 !
 ! NOTE:
 !  If we change size so that it changes with every direction then
 !  we need to change how the code around the put function operates.
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine AB_XFER_finish(xfer, put, put4, ok)
    type (AB_XFER), intent(inout) :: xfer
    logical, intent(out) :: ok
    logical :: tmp_ok
    external put4
    interface
       subroutine put(index,dir,in_array)
         integer, intent(in) :: index,dir
         real, dimension(:),intent(in) :: in_array
       end subroutine put
    end interface
    integer :: ierror
    integer :: i,j,min_buf,max_buf,p,num_rcv,pos
    integer, pointer :: rcv_order(:,:)

    ! initialize status 
    ok=.true.

    ! copy this processors data from snd to rcv buffer
    if (xfer%snd_buf_size(xfer%me) .ne. 0) then
       xfer%rcv_buf(xfer%me)%array=xfer%snd_buf(xfer%me)%array
    endif


    ! wait for all our data to arrive
    call AB_COMM_XCHNG_finish_rcv(xfer%xchng,tmp_ok)
    if (.not. tmp_ok) then
       ok=.false.
       return
    endif

    ! give received data back to user
    do p=0,xfer%max_pn
       pos=1
       num_rcv=xfer%num_nbrs(p)
       rcv_order=>xfer%nbr_order(p)%array
       do i=1,num_rcv
          min_buf=pos
          max_buf=min_buf+xfer%size(rcv_order(3,i))-1
          call put(rcv_order(1,i),rcv_order(2,i), &
                   xfer%rcv_buf(p)%array(min_buf:max_buf))
          pos=max_buf+1
       enddo
    enddo
  end subroutine AB_XFER_finish


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! NAME: AB_XFER_destroy
 !
 ! PURPOSE: destroy a transfer structure.  Freeing buffers, etc.
 !
 ! INPUTS: xfer - transfer structure to be destroyed.
 !
 ! OUTPUTS: ok - status
 !
 ! HISTORY:
 !  10/20/99 Robert Oehmke: created
 !  12/23/99 Robert Oehmke: changed so all communication between
 !                          2 procs is aggeragated into one big message.
 !   5/05/00 Robert Oehmke: allowed msgs in different directions to vary
 !                          in size.
 !   5/11/00 Robert Oehmke: added corner communication.
 !   2/06/01 Robert Oehmke: moved to seperate module
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine AB_XFER_destroy(xfer, ok)
    type (AB_XFER), intent(inout) :: xfer
    logical, intent(out) :: ok
    logical :: tmp_ok
    integer :: ierror
    integer :: p

    ! initialize status 
    ok=.true.

    ! destroy the exchange structure
    call AB_COMM_XCHNG_destroy(xfer%xchng,tmp_ok)
    if (.not. tmp_ok) then
       ok=.false.
       return
    endif

    ! deallocate buffers
    do p=0,xfer%max_pn
       if (xfer%snd_buf_size(p) .ne. 0) then
          deallocate(xfer%snd_buf(p)%array)
       endif

       if (xfer%rcv_buf_size(p) .ne. 0) then
          deallocate(xfer%rcv_buf(p)%array)
       endif
    enddo

    ! deallocate buffer storage
    deallocate(xfer%snd_buf)
    deallocate(xfer%rcv_buf)

    ! deallocate buffer size storage
    deallocate(xfer%snd_buf_size)
    deallocate(xfer%rcv_buf_size)

    ! Deallocate nbr number array
    deallocate(xfer%num_nbrs)

    ! Deallocate nbr order info
    do p=0,xfer%max_pn
       deallocate(xfer%nbr_order(p)%array)
    enddo

    ! Deallocate nbr order info array
    deallocate(xfer%nbr_order)

  end subroutine AB_XFER_destroy

end module AB_XFER_module


