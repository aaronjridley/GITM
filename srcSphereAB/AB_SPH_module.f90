!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!-*- F90 -*- so emacs thinks this is an f90 file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! NAME: AB_SPH_module
!
! PURPOSE: a module which implements the sphere geometry for 
!          AB.
!
! HISTORY:
!  1/05/01 Robert Oehmke: created
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module AB_SPH_module
  use AB_module
  use AB_COMM_module
  use AB_ERROR_module
  implicit none

  ! misc constants
  real, parameter,private  :: pi=3.141592653589793

  type AB_SPH
     private
     integer :: me,np,max_pn
     integer :: num_lat,num_long
     integer :: cpb_long, cpb_lat, cpb_alt,gcn
     type(AB_COMM), pointer :: comm 
     type (AB_GRP),pointer :: grp
     real :: radius, thickness 
     real, dimension(1:3) :: center
     real :: theta_min, theta_max
  end type AB_SPH

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! NAME: AB_SPH_create
 !
 ! PURPOSE: Adds a sphere shaped shell of AB blocks to the passed in grp
 !
 ! INPUTS: grp       - adaptive block group
 !         cpb_long  - the cells per block longitude
 !         cpb_lat   - the cells per block latitude
 !         cpb_alt   - the cells per block altitude
 !         gcn       - the ghostcell number (note gc's don't apply to alt)
 !         num_lat   - number of latitudinal blocks around the sphere
 !         num_long  - number of longitudinal blocks around the sphere.
 !                     For this function, num_long must be one or even.
 !         n_pole_conn - whether or not to connect across the north pole
 !                     TRUE  -> connect 
 !                     FALSE -> otherwise
 !         s_pole_conn - whether or not to connect across the south pole
 !                     TRUE  -> connect 
 !                     FALSE -> otherwise
 !         radius    - radius of the sphere
 !         thickness - thickness of the sphere
 !         center    - center of the sphere
 !         min_lat_ang   - minimum lat angle
 !         max_lat_ang   - maximu lat angle
 !
 ! 
 ! OUTPUTS: sph    - sphere structure
 !          ok     - status code
 !          
 ! HISTORY:
 !  2/02/01
 !
 ! NOTES:
 !  latitude increases as you go from south to north, and
 !  longitude increases as you go west to east. 
 !
 !  For now this just uses the bottom most set of blocks in the group,
 !  eventually this should only use a passed back available set of blocks, so
 !  we can add things to AB_GRPs. 
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine AB_SPH_create(sph,grp, cpb_long,cpb_lat,cpb_alt,gcn,&
                           num_long, num_lat, n_pole_conn, s_pole_conn, &
                           radius, thickness, center, &
                           theta_min, theta_max, min_blk_out, ok)
    type (AB_SPH), intent(out) :: sph
    type (AB_GRP), intent(inout),target :: grp
    integer, intent(in) :: num_lat, num_long,cpb_long,cpb_lat,cpb_alt,gcn
    logical, intent(in) :: n_pole_conn, s_pole_conn
    real, intent(in) :: radius, thickness, theta_min, theta_max
    real, dimension(1:3) :: center
    logical, intent(out) :: ok
    integer, intent(out) :: min_blk_out
    integer :: i,j,lat,long, cross_pole_long
    integer :: num_per_proc,min_blk,max_blk
    integer :: arr_pos,proc,blk, num_blks
    integer :: ierror
    logical :: tmp_ok
    logical :: next_to_north,next_to_south
    integer, dimension(ab_num_nbrs) :: nbrs_p,nbrs_b
    integer :: me,np,max_pn

    ! intialize error returns
    ok=.true.

    ! setup parallel stuff
    sph%comm=>g_comm ! do something better than this when you have time
    call AB_COMM_get_my_proc_num(sph%comm, me)
    call AB_COMM_get_num_procs(sph%comm, np)
    max_pn=np-1

    ! get the number of available blocks from the group
    ! NOTE: I'm assuming that this number is the same across all
    ! processors
    call AB_GRP_get_max_num_blks(grp,num_blks)

    ! error checking
    if ((num_lat*num_long) > (num_blks*np)) then
       ok=.false.
       call AB_ERROR_set("AB_SPH_create","not enough blocks")
       return
    endif

    if ((num_lat*num_long) < np) then
       ok=.false.
       call AB_ERROR_set("AB_SPH_create","too many processors")
       return
    endif

    if ((num_long .ne. 1) .and. (mod(num_long,2) .ne. 0)) then
       ok=.false.
       call AB_ERROR_set("AB_SPH_create","num_long not even")
       return
    endif

    ! Initialize sph variables
    sph%cpb_long=cpb_long
    sph%cpb_lat=cpb_lat
    sph%cpb_alt=cpb_alt
    sph%gcn=gcn
    sph%num_long=num_long
    sph%num_lat=num_lat
    sph%radius=radius
    sph%thickness=thickness
    sph%center=center
    sph%theta_min=theta_min
    sph%theta_max=theta_max
    sph%grp=>grp
    sph%me=me
    sph%np=np
    sph%max_pn=max_pn

    ! NOTE: for the following calcs. lat, long and blks are all numbered
    !       starting with 0 because this makes calc. easier. Outside
    !       this piece of code however they number starting with 1.
    !       This means lat ranges from 0 to num_lat-1, and long ranges
    !       from 0 to num_long-1.

    ! calculate our range of blocks
    call calc_block_range(min_blk,max_blk,num_per_proc,num_lat*num_long, &
                          me,np)

    min_blk_out = min_blk

    ! tell the group how many we're going to use
    call AB_GRP_set_num_used(grp,max_blk-min_blk+1)

    ! Setup new AB blocks
    do i=min_blk,max_blk
       lat=i/num_long
       long=i-(lat*num_long)
       call AB_create(grp,i-min_blk+1,long+1,lat+1)
    enddo

    ! Connect blks into a spherical shell
    do i=min_blk,max_blk
       !! calc position in blk array
       arr_pos=i-min_blk+1

       !! get long and lat of block
       lat=i/num_long
       long=i-(lat*num_long)
     
       !! intialize pole variables 
       next_to_south=.false.
       next_to_north=.false.

       !! do north neighbor connection
       if (lat < num_lat-1) then
          call coord_to_pb(proc,blk,lat+1,long,num_long,num_per_proc,max_pn)
          nbrs_p(ab_north)=proc
          nbrs_b(ab_north)=blk+1
       else
          if (n_pole_conn) then
             if (num_long .ne. 1) then
                cross_pole_long=mod(long+(num_long/2),num_long)
                call coord_to_pb(proc,blk,lat,cross_pole_long, &
                                         num_long, num_per_proc,max_pn)
                nbrs_p(ab_north)=proc
                nbrs_b(ab_north)=blk+1
             else
                nbrs_p(ab_north)=me
                nbrs_b(ab_north)=arr_pos
             endif
          else
             nbrs_p(ab_north)=-1
             nbrs_b(ab_north)=-1
          endif
          next_to_north=.true.
       endif

       !! do south neighbor connection
       if (lat > 0 ) then
          call coord_to_pb(proc,blk,lat-1,long,num_long,num_per_proc,max_pn)
          nbrs_p(ab_south)=proc
          nbrs_b(ab_south)=blk+1
       else
          if (s_pole_conn) then
             if (num_long .ne. 1) then
                cross_pole_long=mod(long+(num_long/2),num_long)
                call coord_to_pb(proc,blk,lat,cross_pole_long, &
                                         num_long, num_per_proc,max_pn)
                nbrs_p(ab_south)=proc
                nbrs_b(ab_south)=blk+1
             else
                nbrs_p(ab_south)=me
                nbrs_b(ab_south)=arr_pos
             endif
          else
             nbrs_p(ab_south)=-1
             nbrs_b(ab_south)=-1
          endif
          next_to_south=.true.
       endif

       !! do east neighbor connection
       if (long < num_long-1 ) then
          call coord_to_pb(proc,blk,lat,long+1,num_long,num_per_proc,max_pn)
          nbrs_p(ab_east)=proc
          nbrs_b(ab_east)=blk+1
       else
          call coord_to_pb(proc,blk,lat,0,num_long,num_per_proc,max_pn)
          nbrs_p(ab_east)=proc
          nbrs_b(ab_east)=blk+1
       endif

       !! do west neighbor connection
       if (long > 0 ) then
          call coord_to_pb(proc,blk,lat,long-1,num_long,num_per_proc,max_pn)
          nbrs_p(ab_west)=proc
          nbrs_b(ab_west)=blk+1
       else
         call coord_to_pb(proc,blk,lat,num_long-1,num_long,num_per_proc,max_pn)
          nbrs_p(ab_west)=proc
          nbrs_b(ab_west)=blk+1
       endif

       !! do northeast corner connection
       if (lat < num_lat-1) then
          if (long < num_long-1 ) then
             call coord_to_pb(proc,blk,lat+1,long+1,num_long, &
                                      num_per_proc,max_pn)
             nbrs_p(ab_northeast)=proc
             nbrs_b(ab_northeast)=blk+1
          else
             call coord_to_pb(proc,blk,lat+1,0,num_long,num_per_proc,max_pn)
             nbrs_p(ab_northeast)=proc
             nbrs_b(ab_northeast)=blk+1
          endif
       else
          if (n_pole_conn) then
             if (num_long .ne. 1) then
                cross_pole_long=mod(long+(num_long/2)+1,num_long)
                call coord_to_pb(proc,blk,lat,cross_pole_long, &
                                         num_long, num_per_proc,max_pn)
                nbrs_p(ab_northeast)=proc
                nbrs_b(ab_northeast)=blk+1
             else
                nbrs_p(ab_northeast)=me
                nbrs_b(ab_northeast)=arr_pos
             endif
          else
             nbrs_p(ab_northeast)=-1
             nbrs_b(ab_northeast)=-1
          endif
       endif

       !! do northwest corner connection
       if (lat < num_lat-1) then
          if (long > 0 ) then
             call coord_to_pb(proc,blk,lat+1,long-1,num_long, &
                                      num_per_proc,max_pn)
             nbrs_p(ab_northwest)=proc
             nbrs_b(ab_northwest)=blk+1
          else
             call coord_to_pb(proc,blk,lat+1,num_long-1,num_long, &
                                      num_per_proc,max_pn)
             nbrs_p(ab_northwest)=proc
             nbrs_b(ab_northwest)=blk+1
          endif
       else
          if (n_pole_conn) then
             if (num_long .ne. 1) then
                ! c_p_l=(long+(num_long/2)-1) mod num_long
                cross_pole_long=mod(long+(num_long/2)+num_long-1,num_long)
                call coord_to_pb(proc,blk,lat,cross_pole_long, &
                                         num_long, num_per_proc,max_pn)
                nbrs_p(ab_northwest)=proc
                nbrs_b(ab_northwest)=blk+1
             else
                nbrs_p(ab_northwest)=me
                nbrs_b(ab_northwest)=arr_pos
             endif
          else
             nbrs_p(ab_northwest)=-1
             nbrs_b(ab_northwest)=-1
          endif
       endif

       !! do southeast corner connection
       if (lat > 0) then
          if (long < num_long-1 ) then
             call coord_to_pb(proc,blk,lat-1,long+1,num_long, &
                                      num_per_proc,max_pn)
             nbrs_p(ab_southeast)=proc
             nbrs_b(ab_southeast)=blk+1
          else
             call coord_to_pb(proc,blk,lat-1,0,num_long,num_per_proc,max_pn)
             nbrs_p(ab_southeast)=proc
             nbrs_b(ab_southeast)=blk+1
          endif
       else
          if (s_pole_conn) then
             if (num_long .ne. 1) then
                cross_pole_long=mod(long+(num_long/2)+1,num_long)
                call coord_to_pb(proc,blk,lat,cross_pole_long, &
                                         num_long, num_per_proc,max_pn)
                nbrs_p(ab_southeast)=proc
                nbrs_b(ab_southeast)=blk+1
             else
                nbrs_p(ab_southeast)=me
                nbrs_b(ab_southeast)=arr_pos
             endif
          else
             nbrs_p(ab_southeast)=-1
             nbrs_b(ab_southeast)=-1
          endif
       endif

       !! do southwest corner connection
       if (lat > 0 ) then
          if (long > 0 ) then
             call coord_to_pb(proc,blk,lat-1,long-1,num_long, &
                                      num_per_proc,max_pn)
             nbrs_p(ab_southwest)=proc
             nbrs_b(ab_southwest)=blk+1
          else
             call coord_to_pb(proc,blk,lat-1,num_long-1,num_long, &
                                      num_per_proc,max_pn)
             nbrs_p(ab_southwest)=proc
             nbrs_b(ab_southwest)=blk+1
          endif
       else
          if (s_pole_conn) then
             if (num_long .ne. 1) then
                ! c_p_l=(long+(num_long/2)-1) mod num_long
                cross_pole_long=mod(long+(num_long/2)+num_long-1,num_long)
                call coord_to_pb(proc,blk,lat,cross_pole_long, &
                                         num_long, num_per_proc,max_pn)
                nbrs_p(ab_southwest)=proc
                nbrs_b(ab_southwest)=blk+1
             else
                nbrs_p(ab_southwest)=me
                nbrs_b(ab_southwest)=arr_pos
             endif
          else
             nbrs_p(ab_southwest)=-1
             nbrs_b(ab_southwest)=-1
          endif
       endif

       ! set the current block to use these neighbor connections
       call AB_set_all_nbrs(grp,arr_pos,nbrs_p,nbrs_b, &
                            next_to_north,next_to_south)
    enddo

    contains
      ! Note: this function assumes ranges start at 0
      !       unlike in the actual structs where the ranges will start at 1
      subroutine calc_block_range(min_blk,max_blk,num_per_proc,num_tot, &
                                  me,np)
        integer, intent(in)  :: num_tot,me,np
        integer, intent(out) :: min_blk,max_blk,num_per_proc

        num_per_proc=num_tot/np

        if (me == np-1) then
           min_blk=num_per_proc*me
           max_blk=num_tot-1
        else
           min_blk=num_per_proc*me
           max_blk=(num_per_proc*(me+1))-1
        endif
        
      end subroutine calc_block_range

      ! Converts long and lat to processor and block number
      ! Note: this function assumes ranges start at 0
      !       unlike in the actual structs where the ranges will start at 1
      subroutine coord_to_pb(proc,blk,lat,long,num_long,num_per_proc, &
                                     max_pn)
        integer, intent(in)  :: lat,long,num_long,num_per_proc,max_pn
        integer, intent(out) :: proc,blk
        integer :: blk_pos
    
        ! convert long and lat into block postion
        blk_pos=(lat*num_long)+long
    
        ! convert the block position into proc and block numbers
        proc=blk_pos/num_per_proc
        if (proc > max_pn) then
           proc=max_pn
        endif
        blk=blk_pos-(proc*num_per_proc)
      end subroutine coord_to_pb


    end subroutine AB_SPH_create

    subroutine AB_SPH_get_xyztp(sph, ind, lat_distr, x, y, z, theta, phi)
      type (AB_SPH), intent(in) :: SPH
      integer, intent(in) :: ind
      real, dimension(1-sph%gcn:,1-sph%gcn:,1:) :: x,y,z,theta,phi
      interface
         subroutine lat_distr(x,y)
           real, intent(in):: x
           real :: y
         end subroutine lat_distr
      end interface
      integer :: long,lat
      integer :: a,i,j
      real :: lgp,th,ph,alt,lout
      integer :: cpb_long,cpb_lat,cpb_alt,gcn
      real :: alt_cell_width
      real :: lat_blk_width,lat_cell_width,lat_gc_width
      real :: phi_blk_width,phi_cell_width,phi_gc_width
      real :: theta_range,theta_min
      real :: cx,cy,cz

      ! get the coordinates of this block
      call AB_get_coords(sph%grp,ind,long,lat)

      ! initialize some handy variables
      cpb_long=sph%cpb_long
      cpb_lat=sph%cpb_lat
      cpb_alt=sph%cpb_alt
      gcn=sph%gcn
      alt_cell_width=sph%thickness/cpb_alt
      phi_blk_width=2*pi/sph%num_long
      phi_cell_width=phi_blk_width/sph%cpb_long
      phi_gc_width=gcn*phi_cell_width
      lat_blk_width=1.0/sph%num_lat
      lat_cell_width=lat_blk_width/cpb_lat
      lat_gc_width=gcn*lat_cell_width
      cx=sph%center(1)
      cy=sph%center(2)
      cz=sph%center(3)
      theta_range=sph%theta_max-sph%theta_min
      theta_min=sph%theta_min

      !! initialize coordinates
      alt=sph%radius
      do a=1,cpb_alt
        lgp=real(lat-1)*lat_blk_width+0.5*lat_cell_width-&
            lat_gc_width        
        do i=1-gcn,cpb_lat+gcn
           ph=real(long-1)*phi_blk_width+0.5*phi_cell_width-phi_gc_width
           do j=1-gcn,cpb_long+gcn

              !!! translate the latitude grid position into theta
              call lat_distr(lgp,lout)
              th=(lout*theta_range)+theta_min

              !!! spherical
              theta(j,i,a)=th
              phi(j,i,a)=ph

              !!! rectangular
              x(j,i,a)=alt*cos(th)*cos(ph+pi)+cx
              y(j,i,a)=alt*cos(th)*sin(ph+pi)+cy
              z(j,i,a)=alt*sin(th)+cz

              ph=ph+phi_cell_width
           enddo
           lgp=lgp+lat_cell_width
        enddo
        alt=alt+alt_cell_width
     enddo

   end subroutine AB_SPH_get_xyztp


end module AB_SPH_module


