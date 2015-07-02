!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine read_waccm_tides

  use ModGitm
  use ModTime
  use ModInputs
  use ModTides
  use ModConstants, only: pi

  implicit none

  logical :: IsFirstTime = .true.
  character(len=40)              :: Waccm_File_Name
  real                           :: version 
  integer                        :: nAlts_gswm, nVars, iVar, iHour
  character(len=40), allocatable :: varName(:)
  integer                        :: iTime(7)
  real, allocatable              :: param(:,:,:,:)
  real, dimension(:,:,:), allocatable :: v24c, v24s, v12c, v12s
  integer :: iError, iBlock, iLon, iLat, iAlt, iLonW, iLatW, i, j
  real :: lon, lat, ri, rj, a
  integer :: iAltW

  call report("read_waccm_tides",1)
  call start_timing("read_waccm_tides")

  if (IsFirstTime) then
     WACCM_file_name ="UA/DataIn/waccm_tides_2000_03.bin" 
     IsFirstTime = .false.
  endif

  iError = 0

  if (iDebugLevel > 2) &
       write(*,*) "===> Reading File : ",WACCM_file_name
  open(iInputUnit_,file=WACCM_file_name,&
       status="old",form="unformatted")

  read(iInputUnit_) version
  read(iInputUnit_) nLonsWaccm, nLatsWaccm, nAltsWaccm
  read(iInputUnit_) nVars
          
  if(.not.allocated(varname)) allocate(varname(nvars))
  
  do ivar = 1, nvars 
     read(iInputUnit_) varname(ivar) 
  enddo
          
  read(iInputUnit_) iTime

  if(.not.allocated(param)) &
       allocate(param(nLonsWaccm, nLatsWaccm, nAltsWaccm, nvars))
        
  if(.not.allocated(t_waccm)) then
     allocate( &
          t_waccm(nLonsWaccm,nLatsWaccm,nAltsWaccm), &
          t0_waccm(nLonsWaccm, nLatsWaccm, nAltsWaccm), &
          ta1_waccm(nLonsWaccm, nLatsWaccm, nAltsWaccm), &
          ta2_waccm(nLonsWaccm, nLatsWaccm, nAltsWaccm), &
          tp1_waccm(nLonsWaccm, nLatsWaccm, nAltsWaccm), &
          tp2_waccm(nLonsWaccm, nLatsWaccm, nAltsWaccm), &
          u_waccm(nLonsWaccm, nLatsWaccm, nAltsWaccm), &
          u0_waccm(nLonsWaccm, nLatsWaccm, nAltsWaccm), &
          ua1_waccm(nLonsWaccm, nLatsWaccm, nAltsWaccm), &
          ua2_waccm(nLonsWaccm, nLatsWaccm, nAltsWaccm), &
          up1_waccm(nLonsWaccm, nLatsWaccm, nAltsWaccm), &
          up2_waccm(nLonsWaccm, nLatsWaccm, nAltsWaccm), &
          v_waccm(nLonsWaccm, nLatsWaccm, nAltsWaccm), &
          v0_waccm(nLonsWaccm, nLatsWaccm, nAltsWaccm), &
          va1_waccm(nLonsWaccm, nLatsWaccm, nAltsWaccm), &
          va2_waccm(nLonsWaccm, nLatsWaccm, nAltsWaccm), &
          vp1_waccm(nLonsWaccm, nLatsWaccm, nAltsWaccm), &
          vp2_waccm(nLonsWaccm, nLatsWaccm, nAltsWaccm), &
          w_waccm(nLonsWaccm, nLatsWaccm, nAltsWaccm), &
          omega_waccm(nLonsWaccm, nLatsWaccm, nAltsWaccm), &
          omega0_waccm(nLonsWaccm, nLatsWaccm, nAltsWaccm), &
          omegaa1_waccm(nLonsWaccm, nLatsWaccm, nAltsWaccm), &
          omegaa2_waccm(nLonsWaccm, nLatsWaccm, nAltsWaccm), &
          omegap1_waccm(nLonsWaccm, nLatsWaccm, nAltsWaccm), &
          omegap2_waccm(nLonsWaccm, nLatsWaccm, nAltsWaccm), &
          v24c(nLonsWaccm, nLatsWaccm, nAltsWaccm), &
          v24s(nLonsWaccm, nLatsWaccm, nAltsWaccm), &
          v12c(nLonsWaccm, nLatsWaccm, nAltsWaccm), &
          v12s(nLonsWaccm, nLatsWaccm, nAltsWaccm), &
          Alt_Waccm_GitmGrid(-1:nLons+2,-1:nLats+2,nAltsWaccm,nBlocks))

     v24c = 0.0
     v24s = 0.0
     v12c = 0.0
     v12s = 0.0
     t_waccm = 0.0
     t0_waccm = 0.0
     ta1_waccm = 0.0
     ta2_waccm = 0.0
     tp1_waccm = 0.0
     tp2_waccm = 0.0
     u_waccm = 0.0
     u0_waccm = 0.0
     ua1_waccm = 0.0
     ua2_waccm = 0.0
     up1_waccm = 0.0
     up2_waccm = 0.0
     v_waccm = 0.0
     v0_waccm = 0.0
     va1_waccm = 0.0
     va2_waccm = 0.0
     vp1_waccm = 0.0
     vp2_waccm = 0.0
     w_waccm = 0.0
     omega_waccm = 0.0
     omega0_waccm = 0.0
     omegaa1_waccm = 0.0
     omegaa2_waccm = 0.0
     omegap1_waccm = 0.0
     omegap2_waccm = 0.0
     Alt_Waccm_GitmGrid = 0.0

     if(.not.allocated(lon_Waccm)) allocate(lon_waccm(nLonsWaccm))
     if(.not.allocated(lat_waccm)) allocate(lat_waccm(nLatsWaccm))
     if(.not.allocated(alt_waccm)) &
          allocate(alt_waccm(nLonsWaccm, nLatsWaccm, nAltsWaccm))
     if(.not.allocated(pressure_waccm)) allocate(pressure_waccm(nAltsWaccm))
        
     do ivar = 1, nvars 
        read(iInputUnit_) param(:,:,:,ivar) 
     enddo

     lon_waccm = param(:,1,1,1)
     lat_waccm = param(1,:,1,2)
     alt_waccm = param(:,:,:,3)
     pressure_waccm = param(1,1,:,4)

     t0_waccm   = param(:,:,:,5)
     v24c = param(:,:,:,6)
     v24s = param(:,:,:,7)
     v12c = param(:,:,:,8)
     v12s = param(:,:,:,9)
     call convert_waccm_fields(v24c, v24s, v12c, v12s,&
          ta1_waccm, ta2_waccm, tp1_waccm, tp2_waccm)

     u0_waccm   = param(:,:,:,10)
     v24c = param(:,:,:,11)
     v24s = param(:,:,:,12)
     v12c = param(:,:,:,13)
     v12s = param(:,:,:,14)
     call convert_waccm_fields(v24c, v24s, v12c, v12s,&
          ua1_waccm, ua2_waccm, up1_waccm, up2_waccm)

     v0_waccm   = param(:,:,:,15)
     v24c = param(:,:,:,16)
     v24s = param(:,:,:,17)
     v12c = param(:,:,:,18)
     v12s = param(:,:,:,19)
     call convert_waccm_fields(v24c, v24s, v12c, v12s,&
          va1_waccm, va2_waccm, vp1_waccm, vp2_waccm)

     omega0_waccm = param(:,:,:,20)
     v24c = param(:,:,:,21)
     v24s = param(:,:,:,22)
     v12c = param(:,:,:,23)
     v12s = param(:,:,:,24)
     call convert_waccm_fields(v24c, v24s, v12c, v12s,&
          omegaa1_waccm, omegaa2_waccm, omegap1_waccm, omegap2_waccm)

  endif

  do iBlock = 1, nBlocks

     do iLon = -1, nLons+2
        lon = mod(longitude(iLon,iBlock)+2*pi,2*pi)
        iLon_Waccm(iLon,iBlock) = 1
        do iLonW = 2, nLonsWaccm-1
           if (lon_waccm(iLonW) <= lon) iLon_Waccm(iLon,iBlock) = iLonW
        enddo
        rLon_Waccm(iLon,iBlock) = &
             (lon - lon_waccm(iLon_Waccm(iLon,iBlock))) / &
             (lon_waccm(iLon_Waccm(iLon,iBlock)+1) - &
              lon_waccm(iLon_Waccm(iLon,iBlock)))
     enddo

     do iLat = -1, nLats+2
        lat = Latitude(iLat,iBlock)
        if (lat >  pi) lat =  2*pi - lat
        if (lat < -pi) lat = -2*pi - lat
        iLat_Waccm(iLat,iBlock) = 1
        do iLatW = 2, nLatsWaccm-1
           if (lat_waccm(iLatW) <= lat) iLat_Waccm(iLat,iBlock) = iLatW
        enddo
        rLat_Waccm(iLat,iBlock) = &
             (lat - lat_waccm(iLat_Waccm(iLat,iBlock))) / &
             (lat_waccm(iLat_Waccm(iLat,iBlock)+1) - &
              lat_waccm(iLat_Waccm(iLat,iBlock)))
     enddo

! This is for testing the linear interpolation
!     do iLon = 1, nLonsWaccm
!        do iLat = 1, nLatsWaccm
!           do iAlt = 1, nAltsWaccm
!              Alt_waccm(iLon,iLat,iAlt) = lat_waccm(iLat)
!           enddo
!        enddo
!     enddo

     do iLon = -1, nLons+2
        do iLat = -1, nLats+2
           i  = iLat_Waccm(iLat,iBlock)
           ri = rLat_Waccm(iLat,iBlock)
           j  = iLon_Waccm(iLon,iBlock)
           rj = rLon_Waccm(iLon,iBlock)
           do iAlt = 1, nAltsWaccm
              ! This interpolates and reverses the altitudes.
              Alt_Waccm_GitmGrid(iLon,iLat,nAltsWaccm-iAlt+1,iBlock) = &
                   (1-rj) * (1-ri) * Alt_waccm(j  ,  i,iAlt) + &
                   (1-rj) * (  ri) * Alt_waccm(j  ,i+1,iAlt) + &
                   (  rj) * (1-ri) * Alt_waccm(j+1,  i,iAlt) + &
                   (  rj) * (  ri) * Alt_waccm(j+1,i+1,iAlt)
           enddo
        enddo
     enddo

     do iLon = -1, nLons+2
        do iLat = -1, nLats+2
           do iAlt = -1, 0
              a = Altitude_gb(iLon,iLat,iAlt,iBlock)
              do iAltW = 2, nAltsWaccm-1
                 if (Alt_Waccm_GitmGrid(iLon,iLat,iAltW,iBlock) <= a) &
                      iAlt_Waccm(iLon,iLat,iAlt+2,iBlock) = iAltW
              enddo
              i = iAlt_Waccm(iLon,iLat,iAlt+2,iBlock)
              rAlt_Waccm(iLon,iLat,iAlt+2,iBlock) = &
                   (a - Alt_waccm_GitmGrid(iLon,iLat,i,iBlock)) / &
                   (Alt_waccm_GitmGrid(iLon,iLat,i+1,iBlock) - &
                    Alt_waccm_GitmGrid(iLon,iLat,i  ,iBlock))
           enddo
        enddo
     enddo

  enddo

  call end_timing("read_waccm_tides")

end subroutine read_waccm_tides

subroutine convert_waccm_fields(v24c, v24s, v12c, v12s, a1, a2, p1, p2)

  use ModTides

  implicit none

  real, dimension(nLonsWaccm, nLatsWaccm, nAltsWaccm), intent(in) :: &
       v24c, v24s, v12c, v12s

  real, dimension(nLonsWaccm, nLatsWaccm, nAltsWaccm), intent(out) :: &
       a1, a2, p1, p2

  a1 = sqrt(v24c**2 + v24s**2)
  a2 = sqrt(v12c**2 + v12s**2)

  p1 = atan2(v24s,v24c)
  p2 = atan2(v12s,v12c)

end subroutine convert_waccm_fields

!-----------------------------------------------------------------------

subroutine update_waccm_tides

  use ModSizeGITM
  use ModTime, only: uTime
  use ModPlanet, only: Rotation_Period
  use ModTides
  use ModConstants, only: pi

  implicit none

  real :: ut1, ut2
  integer :: iLon, iLat, iAlt, i, j, k, iBlock
  real :: ri, rj, rk

  call report("update_waccm_tides",1)
  
  ut1 = utime / Rotation_Period * 2.0 * Pi
  ut2 = utime / Rotation_Period * 4.0 * Pi

  t_waccm = t0_waccm + &
       ta1_waccm * cos(ut1 - tp1_waccm) + &
       ta2_waccm * cos(ut2 - tp2_waccm)

  u_waccm = u0_waccm + &
       ua1_waccm * cos(ut1 - up1_waccm) + &
       ua2_waccm * cos(ut2 - up2_waccm)

  v_waccm = v0_waccm + &
       va1_waccm * cos(ut1 - vp1_waccm) + &
       va2_waccm * cos(ut2 - vp2_waccm)

  omega_waccm = omega0_waccm + &
       omegaa1_waccm * cos(ut1 - omegap1_waccm) + &
       omegaa2_waccm * cos(ut2 - omegap2_waccm)

  do iBlock = 1, nBlocks
     do iLon = -1, nLons+2
        i  = iLon_Waccm(iLon,iBlock)
        ri = rLon_Waccm(iLon,iBlock)
        do iLat = -1, nLats+2
           j  = iLat_Waccm(iLat,iBlock)
           rj = rLat_Waccm(iLat,iBlock)
           do iAlt = -1, 0
              k  = iAlt_Waccm(iLon,iLat,iAlt+2,iBlock)
              rk = rAlt_Waccm(iLon,iLat,iAlt+2,iBlock)

              TidesTemp(iLon,iLat,iAlt+2,iBlock) = &
                   (1-ri)*(1-rj)*(1-rk)*t_waccm(i  ,j  ,k  ) + &
                   (  ri)*(1-rj)*(1-rk)*t_waccm(i+1,j  ,k  ) + &
                   (1-ri)*(  rj)*(1-rk)*t_waccm(i  ,j+1,k  ) + &
                   (  ri)*(  rj)*(1-rk)*t_waccm(i+1,j+1,k  ) + &
                   (1-ri)*(1-rj)*(  rk)*t_waccm(i  ,j  ,k+1) + &
                   (  ri)*(1-rj)*(  rk)*t_waccm(i+1,j  ,k+1) + &
                   (1-ri)*(  rj)*(  rk)*t_waccm(i  ,j+1,k+1) + &
                   (  ri)*(  rj)*(  rk)*t_waccm(i+1,j+1,k+1)

              TidesEast(iLon,iLat,iAlt+2,iBlock) = &
                   (1-ri)*(1-rj)*(1-rk)*u_waccm(i  ,j  ,k  ) + &
                   (  ri)*(1-rj)*(1-rk)*u_waccm(i+1,j  ,k  ) + &
                   (1-ri)*(  rj)*(1-rk)*u_waccm(i  ,j+1,k  ) + &
                   (  ri)*(  rj)*(1-rk)*u_waccm(i+1,j+1,k  ) + &
                   (1-ri)*(1-rj)*(  rk)*u_waccm(i  ,j  ,k+1) + &
                   (  ri)*(1-rj)*(  rk)*u_waccm(i+1,j  ,k+1) + &
                   (1-ri)*(  rj)*(  rk)*u_waccm(i  ,j+1,k+1) + &
                   (  ri)*(  rj)*(  rk)*u_waccm(i+1,j+1,k+1)

              TidesNorth(iLon,iLat,iAlt+2,iBlock) = &
                   (1-ri)*(1-rj)*(1-rk)*v_waccm(i  ,j  ,k  ) + &
                   (  ri)*(1-rj)*(1-rk)*v_waccm(i+1,j  ,k  ) + &
                   (1-ri)*(  rj)*(1-rk)*v_waccm(i  ,j+1,k  ) + &
                   (  ri)*(  rj)*(1-rk)*v_waccm(i+1,j+1,k  ) + &
                   (1-ri)*(1-rj)*(  rk)*v_waccm(i  ,j  ,k+1) + &
                   (  ri)*(1-rj)*(  rk)*v_waccm(i+1,j  ,k+1) + &
                   (1-ri)*(  rj)*(  rk)*v_waccm(i  ,j+1,k+1) + &
                   (  ri)*(  rj)*(  rk)*v_waccm(i+1,j+1,k+1)

              TidesOmega(iLon,iLat,iAlt+2,iBlock) = &
                   (1-ri)*(1-rj)*(1-rk)*omega_waccm(i  ,j  ,k  ) + &
                   (  ri)*(1-rj)*(1-rk)*omega_waccm(i+1,j  ,k  ) + &
                   (1-ri)*(  rj)*(1-rk)*omega_waccm(i  ,j+1,k  ) + &
                   (  ri)*(  rj)*(1-rk)*omega_waccm(i+1,j+1,k  ) + &
                   (1-ri)*(1-rj)*(  rk)*omega_waccm(i  ,j  ,k+1) + &
                   (  ri)*(1-rj)*(  rk)*omega_waccm(i+1,j  ,k+1) + &
                   (1-ri)*(  rj)*(  rk)*omega_waccm(i  ,j+1,k+1) + &
                   (  ri)*(  rj)*(  rk)*omega_waccm(i+1,j+1,k+1)

           enddo
        enddo
     enddo
  enddo

end subroutine update_waccm_tides

subroutine read_tides

  use ModGitm
  use ModTime
  use ModInputs
  use ModTides

  implicit none

  integer :: iComp
  logical :: IsFirstTime = .true.
  real                           :: version 
  integer                        :: nAlts_gswm, nVars, iVar, iHour
  character(len=40), allocatable :: varName(:)
  character(len=2)               :: cMonth
  integer                        :: time(7)
  real, allocatable              :: param(:,:,:,:,:)
  integer :: iError

  call report("read_tides",1)
  call start_timing("read_tides")

  if (IsFirstTime) then
     call i2s(iTimeArray(2), cMonth, 2)
     GSWM_file_name(1) ="UA/DataIn/diur_mig_99km_"//cMonth//".bin" 
     GSWM_file_name(2) ="UA/DataIn/diur_nonmig_99km_"//cMonth//".bin" 
     GSWM_file_name(3) ="UA/DataIn/semidiur_mig_99km_"//cMonth//".bin" 
     GSWM_file_name(4) ="UA/DataIn/semidiur_nonmig_99km_"//cMonth//".bin" 
     IsFirstTime = .false.
  endif

  iError = 0

  do iComp = 1,4

     if (UseGswmComp(iComp)) then

        if (iDebugLevel > 0) &
             write(*,*) "=> Reading File : ",gswm_file_name(iComp)
        open(iInputUnit_,file=gswm_file_name(iComp),&
             status="old",form="unformatted")
          
        read(iInputUnit_) version
        read(iInputUnit_) nLonsGswm, nLatsGswm, nAltsGswm, nvars
          
        if(.not.allocated(varname)) allocate(varname(nvars))
  
        do ivar = 1, nvars 
           read(iInputUnit_) varname(ivar) 
        enddo
          
        read(iInputUnit_) time

        if(.not.allocated(param)) &
             allocate(param(nLonsGswm, nLatsGswm, nAltsGswm, nvars, 24))
        if(.not.allocated(u_gswm)) then
           allocate(u_gswm(nLonsGswm, nLatsGswm, nAltsGswm, 24))
           u_gswm(:,:,:,:) = 0.
        endif
        
        if(.not.allocated(v_gswm)) then
           allocate(v_gswm(nLonsGswm, nLatsGswm, nAltsGswm, 24))
           v_gswm(:,:,:,:) = 0.
        endif

        if(.not.allocated(t_gswm)) then
           allocate(t_gswm(nLonsGswm, nLatsGswm, nAltsGswm, 24))
           t_gswm(:,:,:,:) = 0. 
        endif

        if(.not.allocated(lon_gswm)) allocate(lon_gswm(nLonsGswm))
        if(.not.allocated(lat_gswm)) allocate(lat_gswm(nLatsGswm))
        
        do ihour = 1, 24
           do ivar = 1, nvars 
              read(iInputUnit_) param(:,:,:,ivar,ihour) 
           enddo
        enddo

        lon_gswm(:)     = param(:,1,1,1,1)
        lat_gswm(:)     = param(1,:,1,2,1)
        ! The GSWM grid is uniform, so let's just assume this!
        dLatGswm = lat_gswm(2) - lat_gswm(1)
        dLonGswm = lon_gswm(2) - lon_gswm(1)
        t_gswm(:,:,:,:) = t_gswm(:,:,:,:) + param(:,:,:,4,:) 
        u_gswm(:,:,:,:) = u_gswm(:,:,:,:) + param(:,:,:,5,:) 
        v_gswm(:,:,:,:) = v_gswm(:,:,:,:) + param(:,:,:,6,:) 
        
        close(iInputUnit_)

        deallocate(varname)
        deallocate(param)

     endif

  enddo

  call end_timing("read_tides")

end subroutine read_tides

subroutine update_tides

  use ModGitm
  use ModTime, only: iTimeArray
  use ModInputs, only: iDebugLevel
  use ModTides

  implicit none

  integer :: iLonHalf, iBlock, iLon, iLat, iAlt
  integer :: iLonGSWM, iLatGSWM, iFac1, iFac2
  real :: rfac

  iLonHalf = pi / dLonGswm

  rfac = float(iTimeArray(5))/60.0 + float(iTimeArray(6))/3600.0

  iFac1 = iTimeArray(4)+1            ! this will go from 1-24
  iFac2 = mod(iTimeArray(4)+1,24)+1  ! this will go from 2-25 -> 2-24 back to 1

  do iBlock = 1, nBlocks
     do iAlt = 1,2
        do iLat=-1,nLats+2
           do iLon=-1,nLons+2
           
              iLonGSWM = floor(Longitude(iLon, iBlock) / dLonGSWM)
              iLonGSWM = mod(iLonGSWM + iLonHalf + nLonsGSWM, nLonsGswm) + 1
              iLatGSWM = floor((Latitude(iLat, iBlock) + Pi/2) / dLatGswm) + 1 

              iLatGSWM = max(1,iLatGSWM)
              iLatGSWM = min(nLatsGSWM,iLatGSWM)

              if (iDebugLevel > 5) then
                 write(*,*) "iLon, iLat : ", iLonGSWM, iLatGSWM
                 write(*,*) "lat/lon : ",Latitude(iLat, iBlock)*180/pi, &
                      Longitude(iLon, iBlock)*180/pi
                 write(*,*) "gswm lat/lon : ", &
                      lat_gswm(iLatGSWM)*180/pi, lon_gswm(iLonGSWM)*180/pi
              endif

              TidesTemp(iLon,iLat,iAlt,iBlock) = &
                   (1-rfac) * t_gswm(iLonGSWM, iLatGSWM, 1, iFac1) + &
                   (  rfac) * t_gswm(iLonGSWM, iLatGSWM, 1, iFac2)
              TidesEast(iLon,iLat,iAlt,iBlock) = &
                   (1-rfac) * u_gswm(iLonGSWM, iLatGSWM, 1, iFac1) + &
                   (  rfac) * u_gswm(iLonGSWM, iLatGSWM, 1, iFac2)
              TidesNorth(iLon,iLat,iAlt,iBlock) = &
                   (1-rfac) * v_gswm(iLonGSWM, iLatGSWM, 1, iFac1) + &
                   (  rfac) * v_gswm(iLonGSWM, iLatGSWM, 1, iFac2)

           enddo
        enddo
     enddo
  enddo

end subroutine update_tides
