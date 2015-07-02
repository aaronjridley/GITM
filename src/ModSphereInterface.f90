!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModSphereInterface

  use AB_module
  use AB_COMM_module
  use AB_XFER_module
  use AB_SPH_module
  use AB_ERROR_module
  use AB_XFER_array_util
  use AB_XFER_1blk_util

  use ModGITM
  use ModConstants
  use ModSizeGitm
  use ModPlanet, only: nSpecies
  use ModInputs, only: UseIonAdvection

  implicit none

  integer :: iStartBLK

  ! group structure
  type (AB_GRP) :: uam_grp

  ! Sphere structure
  type (AB_SPH) :: uam_sph

  ! xfer structure
  type (AB_XFER) :: uam_xfer

  ! Sphere description
  !! if there is comm. across poles  
  logical :: n_pole_connect
  logical :: s_pole_connect

  !! handy constants
  real    :: alt_cell_width
  real    :: phi_blk_width
  real    :: phi_cell_width
  real    :: phi_gc_width
  real    :: lat_blk_width
  real    :: lat_cell_width
  real    :: lat_gc_width

  !! total number of cells
  integer :: cells_long
  integer :: cells_lat
  integer :: cells_alt

  !! number of blocks in each direction
  integer            :: blks_long, &
       blks_lat

  real                :: radius, &
       thickness
  real                :: theta_min, &
       theta_max, &
       theta_range

  real, dimension(-1:nLons+2,-1:nLats+2,nAlts) :: x,y,z,theta,phi

  real,dimension(1:3) :: center
  integer, parameter :: gcn=2

  !! cells per block (cpb)
  integer :: cpb_long
  integer :: cpb_lat
  integer :: cpb_alt

  ! UAM type defines
  type UAM_ITER
     private
     type (AB_ITER) :: iter     
  end type UAM_ITER


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! NAME: UAM_module_setup
  !
  ! PURPOSE: setup the UAM module. Must be called before any UAM functions.
  !          This function also allocates and intializes the UAM sphere.
  !
  ! INPUTS: par_context - if parallel=mpi_comm else =whatever
  !         cells_long, cells_lat, cells_alt - dimensions of cell grid
  !         blks_long, blks_lat - blocks to divide the grid into
  !         north_on, south_on  - whether to connect across the poles
  !         blks_lat_ang, min_lat_ang - part of the sphere being represented.
  !         cx,cy,cz            - center of sphere
  !         radius              - radius of sphere
  !         thickness           - thickness of sphere
  !
  !
  ! OUTPUTS: ok - status (optional)
  !
  ! HISTORY:
  !  5/4/00 Robert Oehmke: created
  !
  ! NOTES:
  !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine UAM_module_setup(par_context, &
       icells_long, icells_lat, icells_alt, &
       iblks_long, iblks_lat, &
       imin_lat_ang, imax_lat_ang, &
       inorth_on, isouth_on, &
       icx,icy,icz, &
       iradius, ithickness, &
       ok)

    integer, intent(in) :: par_context
    integer, intent(in) :: icells_long, icells_lat, icells_alt
    integer, intent(in) :: iblks_long, iblks_lat
    real, intent(in)    :: imin_lat_ang, imax_lat_ang
    real, intent(in)    :: icx,icy,icz
    real, intent(in)    :: iradius, ithickness
    logical, intent(in) :: inorth_on, isouth_on
    logical, optional, intent(out) :: ok
    type (AB_ITER) :: iter
    integer :: index,lat,long,s,si,t,i,j,a
    logical :: done, tmp_ok
    real :: lgp,th,ph,alt,lout
    integer :: ierror

    ! intialize error returns
    if (present(ok)) ok=.true.

    ! initialize various parameters
    cells_long=icells_long
    cells_lat=icells_lat
    cells_alt=icells_alt
    blks_long=iblks_long
    blks_lat=iblks_lat
    theta_min=imin_lat_ang
    theta_max=imax_lat_ang
    theta_range=theta_max-theta_min
    center(1)=icx
    center(2)=icy
    center(3)=icz
    radius=iradius
    thickness=ithickness
    n_pole_connect=inorth_on
    s_pole_connect=isouth_on


    !  ! see if we can represent the sphere
    !  !! can we evenly distribute the blocks
    !  if (mod(cells_long,blks_long) .ne. 0) then
    !     if (present(ok)) ok=.false.
    !     call AB_ERROR_set("UAM_module_setup", &
    !         "cells_long not evenly divisable by blks_long")
    !     return
    !  endif
    !  if (mod(cells_lat,blks_lat) .ne. 0) then
    !     if (present(ok)) ok=.false.
    !     call AB_ERROR_set("UAM_module_setup", &
    !         "cells_lat not evenly divisable by blks_lat")
    !     return
    !  endif

    !! can we fit within the blocks
    cpb_long=cells_long  ! /blks_long
    if (cpb_long .gt. nLons) then
       if (present(ok)) ok=.false.
       call AB_ERROR_set("UAM_module_setup", &
            "too many cells per block")
       return
    endif

    cpb_lat=cells_lat  ! /blks_lat
    if (cpb_lat .gt. nLats) then
       if (present(ok)) ok=.false.
       call AB_ERROR_set("UAM_module_setup", &
            "too many cells per block")
       return
    endif

    cpb_alt=cells_alt
    if (cpb_alt .gt. nAlts) then
       if (present(ok)) ok=.false.
       call AB_ERROR_set("UAM_module_setup", &
            "too many cells per block")
       return
    endif

    ! setup various handy constants
    alt_cell_width=thickness/cpb_alt
    phi_blk_width=2*pi/blks_long
    phi_cell_width=phi_blk_width/(cpb_long*blks_long)
    phi_gc_width=gcn*phi_cell_width
    lat_blk_width=1.0/blks_lat
    lat_cell_width=lat_blk_width/cpb_lat
    lat_gc_width=gcn*lat_cell_width

    ! setup AB module
    call AB_module_setup(par_context,tmp_ok)
    if (.not. tmp_ok) then
       if (present(ok)) ok=.false.
       return
    endif

    ! create AB_GRP
    call AB_GRP_create(uam_grp, nBlocksMax, tmp_ok)
    if (.not. tmp_ok) then
       if (present(ok)) ok=.false.
       return
    endif

    ! allocate AB sphere
    call AB_SPH_create(uam_sph,uam_grp,cpb_long,cpb_lat,cpb_alt,gcn, &
         blks_long, blks_lat, n_pole_connect, s_pole_connect, &
         radius, thickness, center, &
         theta_min, theta_max, iStartBLK, tmp_ok)
    if (.not. tmp_ok) then
       if (present(ok)) ok=.false.
       return
    endif


    ! Initialize stuff for sphere
    call AB_ITER_create(iter,uam_grp)
    call AB_ITER_reset(iter,index,done)
    do while (.not. done)

       ! setup coordinates
       call AB_SPH_get_xyztp(uam_sph,index,lat_distr, &
            x(:,:,:), y(:,:,:), z(:,:,:), &
            theta(:,:,:), phi(:,:,:))

       Latitude (:,index) = theta(1,:,1)
       Longitude(:,index) = phi(:,1,1)

       call AB_ITER_next(iter,index,done)

    enddo


  end subroutine UAM_module_setup


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! The following six functions do setup for the ghostcell transfer
!!!! To add a new variable just add a new line in each function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine size_vars_1blk(size)

    integer :: size(ab_num_nbrs), i

    size=0

    ! Rho
    call AB_1blk3_gc_add_size(nLons,nLats,nAlts,2,size)

    ! Velocity
    call AB_1blk4_gc_add_size(nLons,nLats,nAlts,3,2,size)

    ! Species
    call AB_1blk4_gc_add_size(nLons,nLats,nAlts,nSpecies,2,size)

    ! Add Vertical Velocities (nSpecies)
    call AB_1blk4_gc_add_size(nLons,nLats,nAlts,nSpecies,2,size)

    ! Temperature
    call AB_1blk3_gc_add_size(nLons,nLats,nAlts,2,size)

    ! iTemperature
    call AB_1blk3_gc_add_size(nLons,nLats,nAlts,2,size)

    ! eTemperature
    call AB_1blk3_gc_add_size(nLons,nLats,nAlts,2,size)

    ! Ion Species
    if (UseIonAdvection) &
         call AB_1blk4_gc_add_size(nLons,nLats,nAlts,nIons,2,size)

  end subroutine size_vars_1blk


  subroutine pack_vars_1blk(index,dir,pole,out_array)

    integer, intent(in) :: dir,index
    logical, intent(in) :: pole
    real, dimension(:),intent(out) :: out_array
    integer :: p, i

    p=1 ! start packing at position 1 in out_array

    call AB_1blk3_gc_pack(nLons,nLats,nAlts,2, &
         Rho(:,:,1:nAlts,index),dir,pole,p,out_array)

    call AB_1blk4_gc_pack(nLons,nLats,nAlts,3,2, &
         Velocity(:,:,1:nAlts,:,index),dir,pole,p,out_array)

    call AB_1blk4_gc_pack(nLons,nLats,nAlts,nSpecies,2, &
         VerticalVelocity(:,:,1:nAlts,1:nSpecies,index),dir,pole,p,out_array)

    call AB_1blk4_gc_pack(nLons,nLats,nAlts,nSpecies,2, &
         NDensityS(:,:,1:nAlts,1:nSpecies,index),dir,pole,p,out_array)

    call AB_1blk3_gc_pack(nLons,nLats,nAlts,2, &
         Temperature(:,:,1:nAlts,index),dir,pole,p,out_array)

    call AB_1blk3_gc_pack(nLons,nLats,nAlts,2, &
         iTemperature(:,:,1:nAlts,index),dir,pole,p,out_array)

    call AB_1blk3_gc_pack(nLons,nLats,nAlts,2, &
         eTemperature(:,:,1:nAlts,index),dir,pole,p,out_array)

    if (UseIonAdvection) &
         call AB_1blk4_gc_pack(nLons,nLats,nAlts,nIons,2, &
         IDensityS(:,:,1:nAlts,:,index),dir,pole,p,out_array)

  end subroutine pack_vars_1blk

  subroutine unpack_vars_1blk(index,dir,in_array)
    integer, intent(in) :: index,dir
    real, dimension(:), intent(in) :: in_array
    integer :: p,i

    p=1 ! start unpacking at position 1 in in_array

    call AB_1blk3_gc_unpack(nLons,nLats,nAlts,2, &
         Rho(:,:,1:nAlts,index),dir,p,in_array)

    call AB_1blk4_gc_unpack(nLons,nLats,nAlts,3,2, &
         Velocity(:,:,1:nAlts,:,index),dir,p,in_array)

    call AB_1blk4_gc_unpack(nLons,nLats,nAlts,nSpecies,2, &
         VerticalVelocity(:,:,1:nAlts,1:nSpecies,index),dir,p,in_array)

    call AB_1blk4_gc_unpack(nLons,nLats,nAlts,nSpecies,2, &
         NDensityS(:,:,1:nAlts,1:nSpecies,index),dir,p,in_array)

    call AB_1blk3_gc_unpack(nLons,nLats,nAlts,2, &
         Temperature(:,:,1:nAlts,index),dir,p,in_array)

    call AB_1blk3_gc_unpack(nLons,nLats,nAlts,2, &
         iTemperature(:,:,1:nAlts,index),dir,p,in_array)

    call AB_1blk3_gc_unpack(nLons,nLats,nAlts,2, &
         eTemperature(:,:,1:nAlts,index),dir,p,in_array)

    if (UseIonAdvection) &
         call AB_1blk4_gc_unpack(nLons,nLats,nAlts,nIons,2, &
         IDensityS(:,:,1:nAlts,:,index),dir,p,in_array)

  end subroutine unpack_vars_1blk

  subroutine size_vars_nblk(size)
    integer :: size(ab_num_nbrs), i

    size=0
    call AB_array3_gc_add_size(nLons,nLats,nAlts,2,size)

    call AB_array4_gc_add_size(nLons,nLats,nAlts,3,2,size)

    ! Vertical Velocity
    call AB_array4_gc_add_size(nLons,nLats,nAlts,nSpecies,2,size)

    call AB_array4_gc_add_size(nLons,nLats,nAlts,nSpecies,2,size)

    ! temperatures
    call AB_array3_gc_add_size(nLons,nLats,nAlts,2,size)
    call AB_array3_gc_add_size(nLons,nLats,nAlts,2,size)
    call AB_array3_gc_add_size(nLons,nLats,nAlts,2,size)

    if (UseIonAdvection) &
         call AB_array4_gc_add_size(nLons,nLats,nAlts,nIons,2,size)

  end subroutine size_vars_nblk

  subroutine pack_vars_nblk(index,dir,pole,out_array)
    integer, intent(in) :: dir,index
    logical, intent(in) :: pole
    real, dimension(:),intent(out) :: out_array
    integer :: p,i, iIon
    real :: tmpI(-1:nLons+2, -1:nLats+2, 1:nAlts, 1:nIons)
    real :: tmpN(-1:nLons+2, -1:nLats+2, 1:nAlts, 1:nSpecies)

    p=1 ! start packing at position 1 in out_array

    call AB_array3_gc_pack(nLons,nLats,nAlts,2, &
         Rho(:,:,1:nAlts,index),dir,pole,p,out_array)

    call AB_array4_gc_pack(nLons,nLats,nAlts,3,2, &
         Velocity(:,:,1:nAlts,:,index),dir,pole,p,out_array)

    tmpN = VerticalVelocity(:,:,1:nAlts,1:nSpecies,index)
    call AB_array4_gc_pack(nLons,nLats,nAlts,nSpecies,2, &
         tmpN,dir,pole,p,out_array)

    tmpN = NDensityS(:,:,1:nAlts,1:nSpecies,index)
    call AB_array4_gc_pack(nLons,nLats,nAlts,nSpecies,2, &
         tmpN,dir,pole,p,out_array)

    call AB_array3_gc_pack(nLons,nLats,nAlts,2, &
         Temperature(:,:,1:nAlts,index),dir,pole,p,out_array)

    call AB_array3_gc_pack(nLons,nLats,nAlts,2, &
         iTemperature(:,:,1:nAlts,index),dir,pole,p,out_array)

    call AB_array3_gc_pack(nLons,nLats,nAlts,2, &
         eTemperature(:,:,1:nAlts,index),dir,pole,p,out_array)

    if (UseIonAdvection) then
       do iIon = 1, nIons
          tmpI(:,:,:,iIon) = IDensityS(:,:,1:nAlts,iIon,index)
       enddo
       call AB_array4_gc_pack(nLons,nLats,nAlts,nIons,2, &
            tmpI,dir,pole,p,out_array)
    endif
  end subroutine pack_vars_nblk


  subroutine unpack_vars_nblk(index,dir,in_array)
    integer, intent(in) :: index,dir
    real, dimension(:), intent(in) :: in_array
    integer :: p,i, iIon
    real :: tmpI(-1:nLons+2, -1:nLats+2, 1:nAlts, 1:nIons)
    real :: tmpN(-1:nLons+2, -1:nLats+2, 1:nAlts, 1:nSpecies)

    p=1 ! start unpacking at position 1 in in_array

    call AB_array3_gc_unpack(nLons,nLats,nAlts,2, &
         Rho(:,:,1:nAlts,index),dir,p,in_array)

    call AB_array4_gc_unpack(nLons,nLats,nAlts,3,2, &
         Velocity(:,:,1:nAlts,:,index),dir,p,in_array)

    tmpN = VerticalVelocity(:,:,1:nAlts,1:nSpecies,index)
    call AB_array4_gc_unpack(nLons,nLats,nAlts,nSpecies,2, &
         tmpN,dir,p,in_array)
    VerticalVelocity(:,:,1:nAlts,1:nSpecies,index) = tmpN

    tmpN = NDensityS(:,:,1:nAlts,1:nSpecies,index)
    call AB_array4_gc_unpack(nLons,nLats,nAlts,nSpecies,2, &
         tmpN,dir,p,in_array)
    NDensityS(:,:,1:nAlts,1:nSpecies,index) = tmpN

    call AB_array3_gc_unpack(nLons,nLats,nAlts,2, &
         Temperature(:,:,1:nAlts,index),dir,p,in_array)

    call AB_array3_gc_unpack(nLons,nLats,nAlts,2, &
         iTemperature(:,:,1:nAlts,index),dir,p,in_array)

    call AB_array3_gc_unpack(nLons,nLats,nAlts,2, &
         eTemperature(:,:,1:nAlts,index),dir,p,in_array)

    if (UseIonAdvection) then
!       call AB_array4_gc_unpack(nLons,nLats,nAlts,nIons,2, &
!            IDensityS(:,:,1:nAlts,:,index),dir,p,in_array)
       do iIon = 1, nIons
          tmpI(:,:,:,iIon) = IDensityS(:,:,1:nAlts,iIon,index)
       enddo
       call AB_array4_gc_unpack(nLons,nLats,nAlts,nIons,2, &
            tmpI,dir,p,in_array)
       do iIon = 1, nIons
          IDensityS(:,:,1:nAlts,iIon,index) = tmpI(:,:,:,iIon)
       enddo
    endif

  end subroutine unpack_vars_nblk


  ! This function only works where ph is in 0,2*pi and
  ! th is in -pi/2,pi/2, and it assumes both are strictly increasing
  subroutine map_grid(num_ph,num_th,ph,th,from, &
       loc_phi,loc_theta,to)
    integer, intent(in) :: num_ph,num_th
    real, dimension(1:num_ph) :: ph
    real, dimension(1:num_th) :: th
    real, dimension(1:num_ph,1:num_th) :: from
    real, dimension(1-gcn:nLons+gcn, &
         1-gcn:nLats+gcn) :: loc_phi,loc_theta,to
    integer :: i,j
    real :: p,t
    integer :: ind_ph,ind_th
    logical :: fnd_ph, fnd_th
    real :: tmp_to

    ! loop through all points in to grid
    do i=1,cpb_lat
       do j=1,cpb_long

          ! get coords of to grid
          p=loc_phi(j,i) 
          t=loc_theta(j,i) 

          ! Convert coordinates to indices in to grid
          call find_in_list(p, num_ph, ph, ind_ph, fnd_ph)
          call find_in_list(t, num_th, th, ind_th, fnd_th)


          ! if the coordinates are within range then map them
          if (fnd_ph .and. fnd_th) then            
             call interpol(ph(ind_ph),  &
                  th(ind_th), &
                  ph(ind_ph+1),&
                  th(ind_th+1),&
                  from(ind_ph,ind_th), &
                  from(ind_ph,ind_th+1), &
                  from(ind_ph+1,ind_th), &
                  from(ind_ph+1,ind_th+1), &
                  p,t,tmp_to)
             to(j,i)=tmp_to
          endif
       enddo
    enddo

  contains

    ! Returns lower index of range containing the point.
    ! Replace this with a binary search when you have time.
    ! Also, we're going through the input points in a 
    ! specific order, take advantage of that.
    subroutine find_in_list(x, num_list,list, out_ind, fnd)
      real, intent(in) :: x
      integer, intent(in) :: num_list
      real, dimension(1:) :: list
      integer, intent(out) :: out_ind
      logical, intent(out) :: fnd
      integer :: i

      fnd=.false.
      do i=1,num_list-1
         if ((x >= list(i)) .and. (x<list(i+1))) then
            out_ind=i
            fnd=.true.
            return
         endif
      enddo

      ! take care of the case where the point is on the upper edge
      if (x==list(num_list)) then
         out_ind=num_list-1
         fnd=.true.
      endif

    end subroutine find_in_list

  end subroutine map_grid

  subroutine map_reg_grid(min_ph,min_th,d_ph,d_th, &
       num_ph,num_th,from, &
       loc_phi,loc_theta,to)
    real, intent(in) :: min_ph,min_th,d_ph,d_th
    integer, intent(in) :: num_ph,num_th
    real, dimension(1:num_ph,1:num_th) :: from
    real, dimension(1-gcn:nLons+gcn, &
         1-gcn:nLats+gcn) :: loc_phi,loc_theta,to
    integer :: i,j
    real :: ph,th
    integer :: ind_ph,ind_th
    real :: tmp_to

    ! loop through all points in to grid
    do i=1,cpb_lat
       do j=1,cpb_long

          ! get coords of to grid
          ph=loc_phi(j,i) 
          th=loc_theta(j,i) 

          ! Convert coordinates to indices in to grid
          ind_ph=int((ph-min_ph)/d_ph)+1
          ind_th=int((th-min_th)/d_th)+1

          ! if the coordinates are within range then map them
          if (((ind_ph >= 1) .and. (ind_ph <=num_ph)) .and. &
               ((ind_th >= 1) .and. (ind_th <=num_th))) then

             call interpol(min_ph+d_ph*(ind_ph-1), &
                  min_th+d_th*(ind_th-1),&
                  min_ph+d_ph*(ind_ph),&
                  min_th+d_th*(ind_th),&
                  from(ind_ph,ind_th), &
                  from(ind_ph,ind_th+1), &
                  from(ind_ph+1,ind_th), &
                  from(ind_ph+1,ind_th+1), &
                  ph,th,tmp_to)
             to(j,i)=tmp_to
          endif
       enddo
    enddo

  end subroutine map_reg_grid

  ! 0->min, 1->max
  ! assumes x0,x1 and y0,y1 are distinct
  subroutine interpol(x0,y0,x1,y1,q00,q01,q10,q11,x_out,y_out,q_out)
    real, intent(in) :: x0,y0,x1,y1,q00,q01,q10,q11,x_out,y_out
    real, intent(out) :: q_out
    real :: d_x_tot, d_y_tot, d_xy_tot
    real :: d_0x,d_1x
    real :: d_0y,d_1y
    real :: c00,c10,c01,c11

    ! find the x distances
    !! total dist across the cell
    d_x_tot=x1-x0

    !! from the 0 side to the point
    d_0x=x_out-x0

    !! from the 1 side to the point
    d_1x=x1-x_out

    ! find the y distances
    !! total dist across the cell
    d_y_tot=y1-y0

    !! from the 0 side to the point
    d_0y=y_out-y0

    !! from the 1 side to the point
    d_1y=y1-y_out


    ! calc. the product of the total distances
    d_xy_tot=d_x_tot*d_y_tot

    ! calculate the interpolation coeffiecients
    c00=d_1x*d_1y/d_xy_tot
    c10=d_0x*d_1y/d_xy_tot
    c01=d_1x*d_0y/d_xy_tot
    c11=d_0x*d_0y/d_xy_tot

    ! do the interpolation
    q_out=c00*q00+c10*q10+c01*q01+c11*q11

  end subroutine interpol


  ! This is the function used to distribute the range of grid lines
  ! over the range of latitude values. The input value x will be
  ! in the range [0,1] and y, the output value, should also fall 
  ! within this range. Setting this function to y=x will 
  ! distribute the grid lines evenly among the latitude range, less
  ! linear distributions can be achived with other functions. 
  ! 0 is the bottom of both the latitude and grid value ranges, similarly
  ! 1 is the top of both. 

  subroutine lat_distr(x,y)

    real, intent(in):: x
    real :: y, xnorm, xdist, lat
    real, parameter :: pi = 3.141592

    call stretch_grid(x,y)

  end subroutine lat_distr




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! NAME: UAM_ITER_create
  !
  ! PURPOSE: create an iterator structure for the current UAM sphere
  !
  ! INPUTS: 
  ! 
  ! OUTPUTS: iter - the newly created iterator.
  !
  ! HISTORY:
  !  5/5/00 Robert Oehmke: created
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine UAM_ITER_create(r_iter)
    type (UAM_ITER), intent(out) :: r_iter

    call AB_ITER_create(r_iter%iter,uam_grp)

  end subroutine UAM_ITER_create

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! NAME: UAM_ITER_reset
  !
  ! PURPOSE: reset an iter structure to its start
  !
  ! INPUTS: r_iter - iteration structure to be operated upon.
  !
  ! OUTPUTS: index - start index
  !          done  - true if there is nothing to be iterated through, 
  !                  false otherwise. If done=true then index is undefined.
  !
  ! HISTORY:
  !  5/5/00 Robert Oehmke: created
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine UAM_ITER_reset(r_iter,index,done)
    type (UAM_ITER), intent(inout) :: r_iter
    integer, intent(out) :: index
    logical, intent(out) :: done

    call AB_ITER_reset(r_iter%iter,index,done)

  end subroutine UAM_ITER_reset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! NAME: UAM_ITER_next
  !
  ! PURPOSE: advances an iter structure to its next position
  !
  ! INPUTS: r_iter - iteration structure to be operated upon.
  !
  ! OUTPUTS: index - next index
  !          done  - true if there is nothing left to iterate through 
  !                  false otherwise. If done=true then index is undefined.
  !
  ! HISTORY:
  !  5/5/00 Robert Oehmke: created
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine UAM_ITER_next(r_iter,index,done)
    type (UAM_ITER), intent(inout) :: r_iter
    integer, intent(out) :: index
    logical, intent(out) :: done

    call AB_ITER_next(r_iter%iter,index,done)

  end subroutine UAM_ITER_next

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! NAME: UAM_XFER_create
  !
  ! PURPOSE: create a transfer structure. Allocate buffers, initialize 
  !          data structs, etc. 
  !
  ! OUTPUTS: ok - status (optional)
  !
  ! HISTORY:
  !  5/4/00 Robert Oehmke: created
  !
  ! NOTES:
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine UAM_XFER_create(ok)
    logical, optional, intent(out) :: ok
    logical :: tmp_ok
    integer :: size(ab_num_nbrs)


    ! calculate amount being sent in each direction
    if (blks_long==1) then
       call size_vars_1blk(size)
    else
       call size_vars_nblk(size)
    endif

    ! create XFER structure
    call AB_XFER_create(uam_xfer,uam_grp,size,0,tmp_ok)
    if (present(ok)) ok=tmp_ok

  end subroutine UAM_XFER_create

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! NAME: UAM_XFER_start
  !
  ! PURPOSE: start the xfer of ghostcells. Used with UAM finish for
  !          communication hiding.          
  !
  ! OUTPUTS: ok - status (optional) 
  !
  ! HISTORY:
  !  5/4/00 Robert Oehmke: created
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine UAM_XFER_start(ok)
    logical, optional, intent(out) :: ok
    logical :: tmp_ok

    if (blks_long==1) then
       call AB_XFER_start(uam_xfer,pack_vars_1blk,pack_vars_1blk, &
            tmp_ok)
    else
       call AB_XFER_start(uam_xfer,pack_vars_nblk,pack_vars_nblk, &
            tmp_ok)
    endif

    if (present(ok)) ok=tmp_ok

  end subroutine UAM_XFER_start

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! NAME: UAM_XFER_finish
  !
  ! PURPOSE: wait for a real transfer to finish and put the transfered data 
  !          back where the user wants it. Used with UAM start for
  !          communication hiding.          
  !
  ! INPUTS:
  !
  ! OUTPUTS: ok - status (optional) 
  !
  ! HISTORY:
  !  5/4/00 Robert Oehmke: created
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine UAM_XFER_finish(ok)
    logical, optional, intent(out) :: ok
    logical :: tmp_ok

    if (blks_long==1) then
       call AB_XFER_finish(uam_xfer,unpack_vars_1blk,unpack_vars_1blk, &
            tmp_ok)
    else
       call AB_XFER_finish(uam_xfer,unpack_vars_nblk,unpack_vars_nblk, &
            tmp_ok)
    endif

    if (present(ok)) ok=tmp_ok 

  end subroutine UAM_XFER_finish


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! NAME: UAM_XFER_at_once
  !
  ! PURPOSE: Start and finish an UAM xfer. 
  !
  ! INPUTS:
  !
  ! OUTPUTS: ok - status (optional) 
  !
  ! HISTORY:
  !  5/4/00 Robert Oehmke: created
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine UAM_XFER_at_once(ok)
    logical, optional, intent(out) :: ok
    logical tmp_ok

    if (present(ok)) ok=.true.

    if (blks_long==1) then
       call AB_XFER_start(uam_xfer,pack_vars_1blk,pack_vars_1blk, &
            tmp_ok)
    else
       call AB_XFER_start(uam_xfer,pack_vars_nblk,pack_vars_nblk, &
            tmp_ok)
    endif

    if (.not. tmp_ok) then
       if (present(ok)) ok=.false.
       return
    endif

    if (blks_long==1) then
       call AB_XFER_finish(uam_xfer,unpack_vars_1blk,unpack_vars_1blk, &
            tmp_ok)
    else
       call AB_XFER_finish(uam_xfer,unpack_vars_nblk,unpack_vars_nblk, &
            tmp_ok)
    endif

    if (.not. tmp_ok) then
       if (present(ok)) ok=.false.
       return
    endif

  end subroutine UAM_XFER_at_once


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! NAME: UAM_XFER_destroy
  !
  ! PURPOSE: destroy a transfer structure.  Freeing buffers, etc.
  !
  ! INPUTS: 
  !
  ! OUTPUTS:
  !
  ! HISTORY:
  !  5/4/00 Robert Oehmke: created
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine UAM_XFER_destroy(ok)
    logical, optional, intent(out) :: ok
    logical :: tmp_ok

    call AB_XFER_destroy(uam_xfer, tmp_ok)

    if (present(ok)) ok=tmp_ok

  end subroutine UAM_XFER_destroy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! NAME: UAM_map_reg_grid_to_Eta
  !
  ! PURPOSE: Map a regularly spaced grid of data to the UAM grid and write
  !          the interpolated values over the ones already present
  !
  ! INPUTS: min_ph,min_th - spherical coordinates of corner of grid data
  !         d_ph,d_th     - gap between each successive grid line
  !         num_ph,num_th - number of grid lines in each direction
  !         grid(num_ph,num_th)  - grid to be mapped to UAM grid         
  !         alt           - altitude to map to
  !         e             - energy level to map to 
  !
  ! OUTPUTS: ok - status
  !
  ! HISTORY:
  !  02/8/01 Robert Oehmke: created
  !
  ! NOTES:
  !  Need to be aware of fact that polor coords wrap. Add logic to take care
  !  of that.
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine UAM_map_reg_grid_to_Eta(min_ph,min_th,d_ph,d_th, &
       num_ph,num_th,grid,alt,e,ok)
    real, intent(in) :: min_ph,min_th,d_ph,d_th
    integer, intent(in) :: num_ph,num_th
    real, dimension(1:num_ph,1:num_th) :: grid
    integer, intent(in) :: alt,e
    logical, optional, intent(out) :: ok
    type (AB_ITER) :: iter
    logical :: done
    integer :: ind

    ! intialize error returns
    if (present(ok)) ok=.true.

    ! Loop through blocks
    call AB_ITER_create(iter,uam_grp)
    call AB_ITER_reset(iter,ind,done)
    do while (.not. done)

       !! add block level intesection filter here

       !! map grid onto this block
       call map_reg_grid(min_ph,min_th,d_ph,d_th, &
            num_ph,num_th,grid, &
            phi(:,:,alt),theta(:,:,alt), &
            Rho(:,:,alt,ind))

       call AB_ITER_next(iter,ind,done)
    enddo


  end subroutine UAM_map_reg_grid_to_Eta


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! NAME: UAM_map_grid_to_Eta
  !
  ! PURPOSE: Map an irregularly spaced grid of data to the UAM grid and write
  !          the interpolated values over the ones already present
  !
  ! INPUTS: num_ph,num_th - number of grid lines in each direction
  !         ph(num_ph)    - phi coordinates must be in 0,2*pi strictly going up 
  !         th(num_th)    - theta coordinates must be in -pi/2,pi/2 "   "    "
  !         grid(num_ph,num_th)  - grid to be mapped to UAM grid         
  !         alt           - altitude to map to
  !         e             - energy level to map to 
  !
  ! OUTPUTS: ok - status
  !
  ! HISTORY:
  !  02/8/01 Robert Oehmke: created
  !
  ! NOTES:
  !  Need to be aware of fact that polor coords wrap. Add logic to take care
  !  of that. so they can give a grid that crosses 0 longitude
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine UAM_map_grid_to_Eta(num_ph,num_th,ph,th,grid,alt,e,ok)
    integer, intent(in) :: num_ph,num_th
    real, dimension(1:num_ph) :: ph
    real, dimension(1:num_th) :: th
    real, dimension(1:num_ph,1:num_th) :: grid
    integer, intent(in) :: alt,e
    logical, optional, intent(out) :: ok
    type (AB_ITER) :: iter
    logical :: done
    integer :: ind

    ! intialize error returns
    if (present(ok)) ok=.true.

    ! Loop through blocks
    call AB_ITER_create(iter,uam_grp)
    call AB_ITER_reset(iter,ind,done)
    do while (.not. done)

       !! add block level intesection filter here for efficiency's sake


       !! map grid onto this block
       call map_grid(num_ph,num_th,ph,th,grid, &
            phi(:,:,alt),theta(:,:,alt), &
            Rho(:,:,alt,ind))

       call AB_ITER_next(iter,ind,done)
    enddo


  end subroutine UAM_map_grid_to_Eta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! NAME: UAM_write_error
  !
  ! PURPOSE: displays the current error message to the screen. Used if
  !          an UAM function returns false on its ok flag.
  !
  ! INPUTS:
  !
  ! OUTPUTS:
  !
  ! HISTORY:
  !  1/30/01 Robert Oehmke: created
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine UAM_write_error()
    call AB_ERROR_write()
  end subroutine UAM_write_error



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! NAME: UAM_module_takedown
  !
  ! PURPOSE: get rid of the UAM module (deallocate storage, etc) Must be 
  !          called after any UAM functions.
  !
  ! HISTORY:
  !  5/4/00 Robert Oehmke: created
  !
  ! NOTES:
  !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine UAM_module_takedown()

    ! get rid of uam sphere
    call AB_GRP_destroy(uam_grp)

    ! get rid of whole AB shebang
    call AB_module_takedown()

  end subroutine UAM_module_takedown



end module ModSphereInterface
