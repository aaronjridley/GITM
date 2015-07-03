!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!-*- F90 -*- so emacs thinks this is an f90 file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! NAME: AB_module
!
! PURPOSE: a 2D adaptive block structure
!
! HISTORY:
!  10/20/99 Robert Oehmke: created
!
! NOTES:
!   Make the comm, me, np, etc. all contained within the grp struct
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module AB_module
  use AB_COMM_module
  use AB_ARRAY_module
  use AB_ERROR_module
  implicit none

  ! direction number defines
  integer, parameter :: ab_num_cnr_nbrs  = 4, &
                        ab_num_edge_nbrs = 4, &
                        ab_num_nbrs      = 8
                        


  ! direction parmeters
  integer, parameter :: ab_north        = 1, &
                        ab_south        = 2, &
                        ab_east         = 3, &
                        ab_west         = 4, &
                        ab_northeast    = 5, &
                        ab_northwest    = 6, &
                        ab_southeast    = 7, &
                        ab_southwest    = 8


  ! direction indices
  integer, parameter :: ab_north_ind      = 1, &
                        ab_south_ind      = 2, &
                        ab_east_ind       = 3, &
                        ab_west_ind       = 4, &
                        ab_northeast_ind  = 1, &
                        ab_northwest_ind  = 2, &
                        ab_southeast_ind  = 3, &
                        ab_southwest_ind  = 4, &
                        ab_cnr_ind_2_dir  = ab_num_edge_nbrs

  ! refinement level changes
  integer, parameter :: ab_level_ref      = -1, &
                        ab_level_same     =  0, &
                        ab_level_crs      =  1

  ! refined neighbor parameters
  integer, parameter :: ab_rnbr_num    =  2, &
                        ab_rnbr_north  =  1, &
                        ab_rnbr_south  =  2, &
                        ab_rnbr_east   =  1, &
                        ab_rnbr_west   =  2, &
                        ab_rnbr_none   =  1


  ! adaptation parameters
  !! subblock indicators
  integer, private, parameter :: adpt_northeast  = 1, &
                                 adpt_northwest  = 2, &
                                 adpt_southeast  = 3, &
                                 adpt_southwest  = 4
                                   

  !! adpt messages
  integer, private, parameter :: adpt_msg_new_blk  = 1, &
                                 adpt_msg_ref      = 2

  !! adpt message lengths
  integer, private, parameter :: adpt_msg_len_new_blk = &
                                         2*(ab_num_nbrs)+ab_num_edge_nbrs+3, &
                                 adpt_msg_len_ref     = 7 



                        

  ! direction reversal array
  integer, private :: g_rev_dir(ab_num_nbrs) = &
                      (/ab_south,ab_north,ab_west,ab_east, &
                        ab_southwest,ab_southeast,ab_northwest,&
                        ab_northeast/)

  ! corner across pole reversal array
  integer, private :: g_pole_cnr(ab_num_nbrs) = &
                      (/-1,-1,-1,-1, &
                        ab_northwest,ab_northeast,ab_southwest,&
                        ab_southeast/)

  ! direction translation array first index is the pole_type from
  ! AB block the second is the direction.
  integer, private :: g_trans_dir(ab_num_nbrs,0:3)
 
  !! direction translation sub arrays (necessary because you can't
  !! initialize a 2D array
  integer, private :: g_trans_dir0(ab_num_nbrs) = &
   (/ab_south,ab_north,ab_west,ab_east,&
     ab_southwest,ab_southeast,ab_northwest,ab_northeast /)
  integer, private :: g_trans_dir1(ab_num_nbrs) = &
   (/ab_north,ab_north,ab_west,ab_east, &
     ab_northwest,ab_northeast,ab_northwest,ab_northeast/)
  integer, private :: g_trans_dir2(ab_num_nbrs) = &
   (/ab_south,ab_south,ab_west,ab_east, &
     ab_southwest,ab_southeast,ab_southwest,ab_southeast/)
  integer, private :: g_trans_dir3(ab_num_nbrs) = &
   (/ab_north,ab_south,ab_west,ab_east, &
     ab_northwest,ab_northeast,ab_southwest,ab_southeast/)

  ! pole indication array first index is the pole_type from
  ! AB block the second is the direction.
  logical, private :: g_pole(ab_num_nbrs,0:3)
 

  !! pole indication sub arrays (necessary because you can't
  !! initialize a 2D array
  logical, private :: g_pole0(ab_num_nbrs) = &
   (/.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false./)
  logical, private :: g_pole1(ab_num_nbrs) = &
   (/.true.,.false.,.false.,.false.,.true.,.true.,.false.,.false./)
  logical, private :: g_pole2(ab_num_nbrs) = &
   (/.false.,.true.,.false.,.false.,.false.,.false.,.true.,.true./)
  logical, private :: g_pole3(ab_num_nbrs) = &
   (/.true.,.true.,.false.,.false.,.true.,.true.,.true.,.true./)


  ! private global variables
  integer, private :: g_num_procs, g_max_pn, g_me
  integer, private, pointer :: g_buf_pos(:)
  type(AB_COMM), target :: g_comm

  ! type definitions
  type, private :: AB 
     private
     integer, dimension(ab_rnbr_num,ab_num_edge_nbrs) :: nbrs_p
     integer, dimension(ab_rnbr_num,ab_num_edge_nbrs) :: nbrs_b
     integer, dimension(ab_num_edge_nbrs) :: nbrs_l ! nbrs ref. level change
     integer, dimension(ab_num_cnr_nbrs) :: cnrs_p
     integer, dimension(ab_num_cnr_nbrs) :: cnrs_b
     integer :: pole_type
     ! pole_type = 0 -> not at a pole
     ! pole_type = 1 -> next to north pole
     ! pole_type = 2 -> next to south pole
     ! pole_type = 3 -> next to both poles
     integer :: x,y
  end type AB

  type AB_GRP
     private
     integer :: max_num_blks
     integer :: max_used_blks
     integer, pointer :: adpt(:)
     type (AB),pointer :: blks(:)
  end type AB_GRP

  type AB_ITER
     private
     integer :: index
     type (AB_GRP),pointer :: grp     
  end type AB_ITER

  type AB_NBR_ITER
     private
     integer :: index, dir
     type (AB_GRP),pointer :: grp     
  end type AB_NBR_ITER

  ! interfaces
  interface AB_NBR_ITER_reset
     module procedure AB_NBR_ITER_reset_pipd, AB_NBR_ITER_reset_pi
  end interface

  interface AB_NBR_ITER_next
     module procedure AB_NBR_ITER_next_pipd, AB_NBR_ITER_next_pi
  end interface

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! NAME: AB_create
 !
 ! PURPOSE: Initialize an adaptive block for action
 !
 ! INPUTS: grp - the group where the block resides
 !         index  - the index of the block in this group
 !         x - the x coordinate of the specified block (longitude)
 !         y - the y coordinate of the specified block (latitude)
 !
 ! OUTPUTS: 
 !
 ! HISTORY:
 !  2/26/01 Robert Oehmke: created
 !
 ! NOTES: 
 !  In the future I want to move the coordinates from the block to the sph 
 !  module.
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine AB_create(grp,index,x,y)
    type (AB_GRP), intent(inout) :: grp
    integer, intent(in) :: index
    type (AB), pointer :: new_abp
    integer, intent(in) :: x, y
    integer :: i,j

    !get pointer to new adaptive block
    new_abp=>grp%blks(index)

    ! initialize neighbors
    do i=1,ab_num_edge_nbrs
       do j=1,ab_rnbr_num
          new_abp%nbrs_p(j,i)=-1
          new_abp%nbrs_b(j,i)=-1
       enddo
       new_abp%nbrs_l(i)=ab_level_same
    enddo

    ! initialize corners
    do i=1,ab_num_cnr_nbrs
       new_abp%cnrs_p(i)=-1
       new_abp%cnrs_b(i)=-1
    enddo
    
    ! initialize other entries
    new_abp%x=x
    new_abp%y=y
    new_abp%pole_type=0

  end subroutine AB_create

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! NAME: AB_set_nbr
 !
 ! PURPOSE: set the neighbor of an adaptive block in a given direction
 !
 ! INPUTS: grp   - the group where the block resides
 !         index - the index of the block in this group
 !         dir   - the neighbor's direction
 !         pole  - if the neighbor is across a pole
 !         nbr_p - processor of the neighbor
 !         nbr_b - block of the neighbor
 !
 ! OUTPUTS: 
 !
 ! HISTORY:
 !  2/26/01 Robert Oehmke: created
 !
 ! NOTES: 
 !  In the future I want to create another version of this function
 !  that lets the user set refined neighbors also, so multilevel 
 !  structures can be created.
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine AB_set_nbr(grp,index,dir,pole,nbr_p,nbr_b)
    type (AB_GRP), intent(inout) :: grp
    integer, intent(in) :: index,dir,nbr_p,nbr_b
    logical, intent(in) :: pole

    ! set nbr info
    if (dir > ab_num_edge_nbrs) then
       grp%blks(index)%cnrs_p(dir-ab_cnr_ind_2_dir)=nbr_p
       grp%blks(index)%cnrs_b(dir-ab_cnr_ind_2_dir)=nbr_b
    else
       grp%blks(index)%nbrs_p(ab_rnbr_none,dir)=nbr_p
       grp%blks(index)%nbrs_b(ab_rnbr_none,dir)=nbr_b
    endif

    ! update pole type
    !!!!!!! WARNING NEED TO FILL THIS IN BEFORE USING
       
  end subroutine AB_set_nbr


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! NAME: AB_set_all_nbrs
 !
 ! PURPOSE: sets all neighbors of an adaptive block (sometime more 
 !          useful than AB_set_nbr for efficiency reasons.)
 !
 ! INPUTS: grp   - the group where the block resides
 !         index - the index of the block in this group
 !         nbrs_p - processor of the neighbors in a given dir
 !         nbrs_b - block of the neighbor in a give dir
 !         north_pole - are we next to the north pole?
 !         south_pole - are we next to the south pole?
 !
 ! OUTPUTS: 
 !
 ! HISTORY:
 !  2/26/01 Robert Oehmke: created
 !
 ! NOTES: 
 !  In the future I want to create another version of this function
 !  that lets the user set refined neighbors also, so multilevel 
 !  structures can be created.
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine AB_set_all_nbrs(grp,index,nbrs_p,nbrs_b,north_pole,south_pole)
    type (AB_GRP), intent(inout) :: grp
    integer, intent(in) :: index
    integer, dimension(ab_num_nbrs) :: nbrs_p,nbrs_b
    logical, intent(in) :: north_pole, south_pole
    type (AB), pointer :: blk
    integer :: i

    !get pointer to adaptive block
    blk=>grp%blks(index)

    ! set edges
    do i=1,ab_num_edge_nbrs
       blk%nbrs_p(ab_rnbr_none,i)=nbrs_p(i)
       blk%nbrs_b(ab_rnbr_none,i)=nbrs_b(i)
    enddo

    ! set corners
    do i=1,ab_num_cnr_nbrs
       blk%cnrs_p(i)=nbrs_p(i+ab_cnr_ind_2_dir)
       blk%cnrs_b(i)=nbrs_b(i+ab_cnr_ind_2_dir)
    enddo

    ! update pole type
    blk%pole_type=0
    if (north_pole) then
       blk%pole_type=blk%pole_type+1
    endif
    if (south_pole) then
       blk%pole_type=blk%pole_type+2
    endif

  end subroutine AB_set_all_nbrs


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! NAME: AB_module_setup
 !
 ! PURPOSE: setup the AB module. Must be called before any AB functions
 !
 ! OUTPUTS: ok - status code
 !
 ! HISTORY:
 !  10/20/99 Robert Oehmke: created
 !  12/23/99 Robert Oehmke: added allocation of global arrays
 !
 ! NOTES:
 !  Perhaps we should get rid of this function and just call the MPI 
 !  initialization functions instead. Actually, it'd be nice to create
 !  a seperate comm for this module so we don't have to worry about 
 !  tag collisions. 
 !  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine AB_module_setup(context, ok)
    integer, intent(in) :: context
    logical, intent(out) :: ok
    logical :: tmp_ok
    integer :: ierror

    ! intialize error returns
    ok=.true.

    ! create comm for this group to use
    call AB_COMM_create(g_comm,context)

    ! set variables
    call AB_COMM_get_my_proc_num(g_comm, g_me)
    call AB_COMM_get_num_procs(g_comm, g_num_procs)
    g_max_pn=g_num_procs-1

    ! initialize direction translation array
    g_trans_dir(:,0)=g_trans_dir0(:)
    g_trans_dir(:,1)=g_trans_dir1(:)
    g_trans_dir(:,2)=g_trans_dir2(:)
    g_trans_dir(:,3)=g_trans_dir3(:)

    ! initialize pole indication array
    g_pole(:,0)=g_pole0(:)
    g_pole(:,1)=g_pole1(:)
    g_pole(:,2)=g_pole2(:)
    g_pole(:,3)=g_pole3(:)


    ! allocate global g_buf_pos array used in transfer procedures
    allocate(g_buf_pos(0:g_max_pn), stat=ierror)
    if (ierror .ne. 0) then
       ok=.false.
       call AB_ERROR_set("AB_module_setup","allocate error ",ierror)
       return
    endif

  end subroutine AB_module_setup


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! NAME: AB_GRP_create_sphere
 !
 ! PURPOSE: create a group of AB blocks
 !
 ! INPUTS: num_blks  - maximum number of blocks that this group can hold.
 !
 ! 
 ! OUTPUTS: grp    - newly created group of adaptive blocks
 !          ok     - status code
 !          
 ! HISTORY:
 !  2/06/01 created: Robert Oehmke
 !
 ! NOTES:
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine AB_GRP_create(grp, num_blks, ok)
    type (AB_GRP), intent(out) :: grp
    integer, intent(in) :: num_blks
    logical, intent(out) :: ok
    integer :: ierror
    logical :: tmp_ok

    ! intialize error returns
    ok=.true.

    ! Allocate space for the blocks
    allocate(grp%blks(num_blks), stat=ierror)
    if (ierror .ne. 0) then
       ok=.false.
       call AB_ERROR_set("AB_GRP_create","allocate error ",ierror)
       return
    endif

    ! Allocate space for adpt values
    allocate(grp%adpt(num_blks), stat=ierror)
    if (ierror .ne. 0) then
       ok=.false.
       call AB_ERROR_set("AB_GRP_create_sphere","allocate error ",ierror)
       return
    endif

    ! Initialize adaptation values
    grp%adpt=0

    ! Initialize variables
    grp%max_num_blks=num_blks   
    grp%max_used_blks=0

  end subroutine AB_GRP_create

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! NAME: AB_GRP_get_max_num_blks
 !
 ! PURPOSE: Get the maximum number of blocks in this group
 !
 ! INPUTS: grp - the ab group of blocks. 
 ! 
 ! OUTPUTS: max_num - the maximum number of blocks in this group
 !
 ! HISTORY:
 !  2/26/01 Robert Oehmke: created
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine AB_GRP_get_max_num_blks(grp,max_num)
    type (AB_GRP), intent(in),target :: grp
    integer, intent(out) :: max_num

    max_num=grp%max_num_blks

  end subroutine AB_GRP_get_max_num_blks

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! NAME: AB_GRP_set_num_used
 !
 ! PURPOSE: Set the number of used blocks in this group 
 !
 ! INPUTS: grp      - the ab group of blocks. 
 !         num_used - the number of used blocks in this group 
 !
 ! OUTPUTS: 
 !
 ! HISTORY:
 !  2/26/01 Robert Oehmke: created
 !
 ! NOTES: 
 !  This function will eventually go away once I have a more
 !  rational system for obtaining blocks from a group.
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine AB_GRP_set_num_used(grp,num_used)
    type (AB_GRP), intent(inout),target :: grp
    integer, intent(in) :: num_used

    grp%max_used_blks=num_used

  end subroutine AB_GRP_set_num_used

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! NAME: AB_ITER_create
 !
 ! PURPOSE: create an iterator structure for the AB group grp
 !
 ! INPUTS: grp - the ab group of blocks over which this iterator will
 !                  range.
 ! 
 ! OUTPUTS: iter - the newly created iterator.
 !
 ! HISTORY:
 !  10/20/99 Robert Oehmke: created
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine AB_ITER_create(iter,grp)
    type (AB_GRP), intent(in),target :: grp
    type (AB_ITER), intent(out) :: iter

    iter%grp=>grp
    iter%index=1

  end subroutine AB_ITER_create

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! NAME: AB_ITER_reset
 !
 ! PURPOSE: reset an iter structure to its start
 !
 ! INPUTS: iter - iteration structure to be operated upon.
 !
 ! OUTPUTS: index - start index
 !          done  - true if there is nothing to be iterated through, 
 !                  false otherwise. If done=true then index is undefined.
 !
 ! HISTORY:
 !  10/20/99 Robert Oehmke: created
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine AB_ITER_reset(iter,index,done)
    type (AB_ITER), intent(inout) :: iter
    integer, intent(out) :: index
    logical, intent(out) :: done

    iter%index=1
    
    if (iter%grp%max_used_blks==0) then
       done=.true.
       index=0
    else
       done=.false.
       index=1
    endif
  end subroutine AB_ITER_reset

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! NAME: AB_ITER_next
 !
 ! PURPOSE: advances an iter structure to its next position
 !
 ! INPUTS: iter - iteration structure to be operated upon.
 !
 ! OUTPUTS: index - next index
 !          done  - true if there is nothing left to iterate through 
 !                  false otherwise. If done=true then index is undefined.
 !
 ! HISTORY:
 !  10/20/99 Robert Oehmke: created
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine AB_ITER_next(iter,index,done)
    type (AB_ITER), intent(inout) :: iter
    integer, intent(out) :: index
    logical, intent(out) :: done

    iter%index=iter%index+1

    if (iter%index>iter%grp%max_used_blks) then
       done=.true.
       index=0
    else
       done=.false.
       index=iter%index
    endif

  end subroutine AB_ITER_next


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! NAME: AB_NBR_ITER_create
 !
 ! PURPOSE: create a neighbor iterator structure for the AB group grp.
 !          A neighbor iterator allows the user to iterate through
 !          all the blocks neighbors. This iteration is gaurenteed
 !          to occur in the same order as the last one, unless
 !          a modification occurs. Modifications can be detected
 !          with the modification signaler:
 !
 ! INPUTS: grp - the ab group of blocks over which this iterator will
 !                  range.
 ! 
 ! OUTPUTS: iter - the newly created neighbor iterator.
 !
 ! HISTORY:
 !  2/2/01 Robert Oehmke: created
 !
 ! NOTES:
 !  Should I only iterate through valid neighbors (e.g. no -1's)?
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine AB_NBR_ITER_create(iter,grp)
    type (AB_GRP), intent(in),target :: grp
    type (AB_NBR_ITER), intent(out) :: iter

    iter%grp=>grp
    iter%index=1
    iter%dir=1

  end subroutine AB_NBR_ITER_create

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! NAME: AB_NBR_ITER_reset_pipd
 !
 ! PURPOSE: reset a neigbor iter structure to its start
 !
 ! INPUTS: iter - iteration structure to be operated upon.
 !
 ! OUTPUTS: index - start index
 !          dir   - start dir
 !          nbr_proc  -  neigbor's processor (-1 if no neighbor)
 !          nbr_index -  neigbor's index  (-1 if no neighbor)
 !          pole   - if our neighbor is across a pole
 !          nbr_dir - the direction we are from our neighbor    
 !          done  - true if there is nothing to be iterated through, 
 !                  false otherwise. If done=true then other outputs are
 !                  undefined.
 !
 ! HISTORY:
 !  02/02/01 Robert Oehmke: created
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine AB_NBR_ITER_reset_pipd(iter,index,dir, &
                               nbr_proc,nbr_index,pole,nbr_dir, &
                               done)
    type (AB_NBR_ITER), intent(inout) :: iter
    integer, intent(out) :: index,dir,nbr_proc,nbr_index,nbr_dir
    logical, intent(out) :: done,pole
    integer :: pole_type

    iter%index=1
    iter%dir=1
    
    if (iter%grp%max_used_blks==0) then
       done=.true.
       index=0
       dir=0
       nbr_proc=-1
       nbr_index=0
       nbr_dir=0
       pole=.false.
    else
       done=.false.
       index=1
       dir=1
       nbr_proc=iter%grp%blks(index)%nbrs_p(ab_rnbr_none,dir)
       nbr_index=iter%grp%blks(index)%nbrs_b(ab_rnbr_none,dir)
       pole_type=iter%grp%blks(index)%pole_type
       pole=g_pole(dir,pole_type)
       nbr_dir=g_trans_dir(dir,pole_type)
    endif
          
  end subroutine AB_NBR_ITER_reset_pipd

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! NAME: AB_NBR_ITER_next_pipd
 !
 ! PURPOSE: advances a neighbor iter structure to its next position
 !
 ! INPUTS: iter - iteration structure to be operated upon.
 !
 ! OUTPUTS: index - start index
 !          dir   - start dir
 !          nbr_proc  -  neigbor's processor (-1 if no neighbor)
 !          nbr_index -  neigbor's index  (-1 if no neighbor)
 !          pole   - if our neighbor is across a pole
 !          nbr_dir - the direction we are from our neighbor     
 !          done  - true if there is nothing to be iterated through, 
 !                  false otherwise. If done=true then other outputs are
 !                  undefined.
 !
 ! HISTORY:
 !  02/02/01 Robert Oehmke: created
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine AB_NBR_ITER_next_pipd(iter,index,dir, &
                                   nbr_proc,nbr_index,pole,nbr_dir, &
                                   done)
    type (AB_NBR_ITER), intent(inout) :: iter
    integer, intent(out) :: index,dir,nbr_proc,nbr_index,nbr_dir
    logical, intent(out) :: done,pole
    integer :: pole_type

    ! advance indicies
    iter%dir=iter%dir+1    
    if (iter%dir>ab_num_nbrs) then
       iter%dir=1
       iter%index=iter%index+1
    endif


    ! check if we're still in range and if so then do output
    if (iter%index>iter%grp%max_used_blks) then
       done=.true.
       index=0
       dir=0
       nbr_proc=-1
       nbr_index=0
       nbr_dir=0
       pole=.false.
    else
       done=.false.
       index=iter%index
       dir=iter%dir
       if (dir > ab_num_edge_nbrs) then
          nbr_proc=iter%grp%blks(index)%cnrs_p(dir-ab_cnr_ind_2_dir)
          nbr_index=iter%grp%blks(index)%cnrs_b(dir-ab_cnr_ind_2_dir)
       else
          nbr_proc=iter%grp%blks(index)%nbrs_p(ab_rnbr_none,dir)
          nbr_index=iter%grp%blks(index)%nbrs_b(ab_rnbr_none,dir)
       endif
       pole_type=iter%grp%blks(index)%pole_type
       pole=g_pole(dir,pole_type)
       nbr_dir=g_trans_dir(dir,pole_type)
    endif
  end subroutine AB_NBR_ITER_next_pipd


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! NAME: AB_NBR_ITER_reset_pi
 !
 ! PURPOSE: reset a neigbor iter structure to its start
 !
 ! INPUTS: iter - iteration structure to be operated upon.
 !
 ! OUTPUTS: index - start index
 !          dir   - start dir
 !          nbr_proc  -  neigbor's processor (-1 if no neighbor)
 !          nbr_index -  neigbor's index  (-1 if no neighbor)
 !          done  - true if there is nothing to be iterated through, 
 !                  false otherwise. If done=true then other outputs are
 !                  undefined.
 !
 ! HISTORY:
 !  02/02/01 Robert Oehmke: created
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine AB_NBR_ITER_reset_pi(iter,index,dir, &
                               nbr_proc,nbr_index, &
                               done)
    type (AB_NBR_ITER), intent(inout) :: iter
    integer, intent(out) :: index,dir,nbr_proc,nbr_index
    logical, intent(out) :: done

    iter%index=1
    iter%dir=1
    
    if (iter%grp%max_used_blks==0) then
       done=.true.
       index=0
       dir=0
       nbr_proc=-1
       nbr_index=0
    else
       done=.false.
       index=1
       dir=1
       nbr_proc=iter%grp%blks(index)%nbrs_p(ab_rnbr_none,dir)
       nbr_index=iter%grp%blks(index)%nbrs_b(ab_rnbr_none,dir)
    endif
          
  end subroutine AB_NBR_ITER_reset_pi

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! NAME: AB_NBR_ITER_next_pi
 !
 ! PURPOSE: advances a neighbor iter structure to its next position
 !
 ! INPUTS: iter - iteration structure to be operated upon.
 !
 ! OUTPUTS: index - start index
 !          dir   - start dir
 !          nbr_proc  -  neigbor's processor (-1 if no neighbor)
 !          nbr_index -  neigbor's index  (-1 if no neighbor)
 !          done  - true if there is nothing to be iterated through, 
 !                  false otherwise. If done=true then other outputs are
 !                  undefined.
 !
 ! HISTORY:
 !  02/02/01 Robert Oehmke: created
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine AB_NBR_ITER_next_pi(iter,index,dir, &
                               nbr_proc,nbr_index, &
                               done)
    type (AB_NBR_ITER), intent(inout) :: iter
    integer, intent(out) :: index,dir,nbr_proc,nbr_index
    logical, intent(out) :: done

    ! advance indicies
    iter%dir=iter%dir+1    
    if (iter%dir>ab_num_nbrs) then
       iter%dir=1
       iter%index=iter%index+1
    endif


    ! check if we're still in range and if so then do output
    if (iter%index>iter%grp%max_used_blks) then
       done=.true.
       index=0
       dir=0
       nbr_proc=-1
       nbr_index=0
    else
       done=.false.
       index=iter%index
       dir=iter%dir
       if (dir > ab_num_edge_nbrs) then
          nbr_proc=iter%grp%blks(index)%cnrs_p(dir-ab_cnr_ind_2_dir)
          nbr_index=iter%grp%blks(index)%cnrs_b(dir-ab_cnr_ind_2_dir)
       else
          nbr_proc=iter%grp%blks(index)%nbrs_p(ab_rnbr_none,dir)
          nbr_index=iter%grp%blks(index)%nbrs_b(ab_rnbr_none,dir)
       endif
    endif
  end subroutine AB_NBR_ITER_next_pi

 subroutine test_NBR_ITER(grp)
   type (AB_GRP), intent(in) :: grp
   type (AB_NBR_ITER) :: iter
   integer :: b,d,n_p,n_b, i
   logical :: done

   i=1
   write(*,*) "AB_NBR_ITER"
   call AB_NBR_ITER_create(iter,grp)   
   call AB_NBR_ITER_reset(iter,b,d,n_p,n_b,done)   
   do while(.not. done)

      write(*,*) i,"[",b,",",d,"]={",n_p,",",n_b,"}"
      i=i+1
      call AB_NBR_ITER_next(iter,b,d,n_p,n_b,done)   
   enddo

   i=1
   write(*,*) 
   write(*,*) "Just Loops" 
   do b=1,grp%max_used_blks

       !! do direct neighbors
       do d=1,ab_num_edge_nbrs
          n_p=grp%blks(b)%nbrs_p(ab_rnbr_none,d)
          n_b=grp%blks(b)%nbrs_b(ab_rnbr_none,d)
          write(*,*) i,"[",b,",",d,"]={",n_p,",",n_b,"}"          
          i=i+1
       enddo

       !! do corner neighbors
       do d=1,ab_num_cnr_nbrs
          n_p=grp%blks(b)%cnrs_p(d)
          n_b=grp%blks(b)%cnrs_b(d)
          write(*,*) i,"[",b,",",d+4,"]={",n_p,",",n_b,"}"
          i=i+1
       enddo
    enddo    




 end subroutine test_NBR_ITER


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! NAME: AB_get_coords
 !
 ! PURPOSE: given an index number returns the coords associated with
 !          that index
 !
 ! INPUTS: grp - the group from which to return the indexed blocks
 !                  coordinates.
 !         index  - the index of the block in this group whose coordinates
 !                  are being requested.
 !
 ! OUTPUTS: x - the x coordinate of the specified block (longitude)
 !          y - the y coordinate of the specified block (latitude)
 !
 ! HISTORY:
 !  10/20/99 Robert Oehmke: created
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine AB_get_coords(grp,index,x,y)
    type (AB_GRP), intent(in) :: grp
    integer, intent(in) :: index
    integer, intent(out) :: x,y
    integer :: f

    x=grp%blks(index)%x
    y=grp%blks(index)%y

  end subroutine AB_get_coords


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! NAME: AB_GRP_destroy
 !
 ! PURPOSE: do all the cleanup neccessary to stop using a AB_GRP structure
 !          (e.g. deallocate all associated memory)
 !
 ! INPUTS: grp - group of adaptive blocks to be destroyed.
 !
 ! OUTPUTS:
 !
 ! HISTORY:
 !  10/20/99 Robert Oehmke: created
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine AB_GRP_destroy(grp)
    type (AB_GRP), intent(inout) :: grp
    integer :: p

    ! Deallocate adpt array
    deallocate(grp%adpt)

    ! Deallocate blk array
    deallocate(grp%blks)

  end subroutine AB_GRP_destroy

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! NAME: AB_get_nbrs
 !
 ! PURPOSE: Returns a list of the neighbors (including corners) of a 
 !          block
 !
 ! INPUTS: grp - group of adaptive blocks.
 !         index  - index of the block to return the neighbors of
 !
 ! OUTPUTS: nbrs - 2xab_num_nbrs array where the first index is 1 for blk #
 !                 and 2 for proc #. and the second index indicates the 
 !                 direction.
 !
 ! HISTORY:
 !  5/11/00 Robert Oehmke: created
 !  
 ! NOTES:
 !  This is mostly for debugging.
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine AB_get_nbrs(grp,index,nbrs)
    type (AB_GRP), intent(in) :: grp
    integer, intent(in) :: index
    integer, intent(out) :: nbrs(2,ab_num_nbrs)
    integer :: i

    do i=1,ab_num_edge_nbrs
       nbrs(1,i)=grp%blks(index)%nbrs_b(ab_rnbr_none,i)
       nbrs(2,i)=grp%blks(index)%nbrs_p(ab_rnbr_none,i)
    enddo

    do i=1,ab_num_cnr_nbrs
       nbrs(1,i+ab_cnr_ind_2_dir)=grp%blks(index)%cnrs_b(i)
       nbrs(2,i+ab_cnr_ind_2_dir)=grp%blks(index)%cnrs_p(i)
    enddo

  end subroutine AB_get_nbrs


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! NAME: AB_module_takedown
 !
 ! PURPOSE: get rid of the AB module (deallocate storage, etc) Must be 
 !          called after any AB functions.
 !
 ! HISTORY:
 !  12/23/99 Robert Oehmke: created
 !
 ! NOTES:
 !  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine AB_module_takedown()
    integer :: ierror
    logical :: tmp_ok

    ! deallocate global g_buf_pos array
    deallocate(g_buf_pos)

  end subroutine AB_module_takedown


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! NAME: AB_set_adpt_value()
 !
 ! PURPOSE: set an adaptation value for a block
 !
 ! INPUTS: grp - the ab group of blocks over in which to set the adpt value.
 !         ind - the block to set the value for
 !         a   - the adapt value
 !
 ! OUTPUTS: none
 !
 ! HISTORY:
 !  2/06/01 Robert Oehmke: created
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine AB_set_adpt_value(grp,ind,a)
    type (AB_GRP), intent(inout) :: grp
    integer, intent(in) :: ind,a

    grp%adpt(ind)=a

  end subroutine AB_set_adpt_value

end module AB_module
