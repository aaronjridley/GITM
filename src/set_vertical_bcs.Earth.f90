!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!\
! ------------------------------------------------------------
! set_boundary
! ------------------------------------------------------------
!/

subroutine set_vertical_bcs(LogRho,LogNS,Vel_GD,Temp, LogINS, iVel, VertVel)

  ! Fill in ghost cells at the top and bottom

  use ModSizeGitm, only: nAlts
  use ModPlanet, only: nSpecies, nIonsAdvect, Mass, nIons, IsEarth,iN2_,iNO_
  use ModGITM, only: TempUnit, iEast_, iNorth_, iUp_
  use ModInputs
  use ModConstants
  use ModTime, only: UTime, iJulianDay,currenttime
  use ModVertical, only: &
       Lat, Lon, Gravity_G, Altitude_G, dAlt_F, iLon1D, iLat1D, iBlock1D, SZAVertical
  use ModIndicesInterfaces, only: get_HPI
  use ModTides, only: TidesNorth, TidesEast, TidesTemp

  use EUA_ModMsis90, ONLY: meter6

  implicit none

  real, intent(inout) :: &
       LogRho(-1:nAlts+2), &
       LogNS(-1:nAlts+2,nSpecies), &
       LogINS(-1:nAlts+2,nIons), &
       Vel_GD(-1:nAlts+2,3), &
       IVel(-1:nAlts+2,3), &
       Temp(-1:nAlts+2), &
       VertVel(-1:nAlts+2,nSpecies)

  integer :: iSpecies, iAlt
  real    :: InvScaleHeightS, InvScaleHgt, Alt, Lst, Ap = 4.0, dn, dt
  logical :: IsFirstTime = .true., UseMsisBCs = .false.
  real    :: HP, v(2)
  integer :: ierror
  real    :: temptemp
  real    :: fac

  integer, dimension(25) :: sw

  !-----------------------------------------------------------
  ! Bottom
  !-----------------------------------------------------------

  if (IsEarth) UseMsisBCs = UseMsis

  if (IsFirstTime .and. UseMsisBCs) then
     call meter6(.true.)
     sw = 1
     IsFirstTime = .true.
  endif

  if (UseMsisBCs) then
     call get_HPI(CurrentTime, HP, iError)  
     if (iError > 0) hp = 40.0
     Ap = min(200.,max(-40.72 + 1.3 * HP, 10.))
     do iAlt = -1, 0
        Alt = Altitude_G(iAlt)/1000.0
        Lst = mod(UTime/3600.0+Lon/15.0,24.0)

        call msis_bcs(iJulianDay,UTime,Alt,Lat,Lon,Lst, &
             F107A,F107,AP,LogNS(iAlt,:), temptemp, &
             LogRho(iAlt),v)

        if (.not. DuringPerturb) temp(iAlt) = temptemp

        vel_gd(iAlt,iEast_) = v(iEast_)
        vel_gd(iAlt,iNorth_) = v(iNorth_)

     enddo
  else
     ! Don't Let the winds blow
     Vel_GD(-1:0,iEast_)  = 0.0
     Vel_GD(-1:0,iNorth_) = 0.0
  endif

  if (.not. DuringPerturb) then
     Vel_GD(-1:0,iUp_)    = 0.0
     VertVel(-1:0,:)      = 0.0
  endif

  if (UseGSWMTides) then
     Vel_GD(-1:0,iEast_)  = TidesEast(iLon1D,iLat1D,1:2,iBlock1D)
     Vel_GD(-1:0,iNorth_) = TidesNorth(iLon1D,iLat1D,1:2,iBlock1D)
     Temp(-1:0)           = TidesTemp(iLon1D,iLat1D,1:2,iBlock1D) + Temp(-1:0)
  endif
  if (UseWACCMTides) then
     Vel_GD(-1:0,iEast_)  = TidesEast(iLon1D,iLat1D,1:2,iBlock1D)
     Vel_GD(-1:0,iNorth_) = TidesNorth(iLon1D,iLat1D,1:2,iBlock1D)
     Temp(-1:0)           = TidesTemp(iLon1D,iLat1D,1:2,iBlock1D)
  endif

  IVel(-1:0,iUp_)      = 0.0

  do iSpecies=1,nIons-1 !Advect
     dn = (LogINS(2,iSpecies) - LogINS(1,iSpecies))
     LogINS(0,iSpecies) = LogINS(1,iSpecies) - dn
     LogINS(-1,iSpecies) = LogINS(0,iSpecies) - dn
  enddo

!  ! Lower boundary for NO on Earth
!  if (nSpecies == iNO_) then
!     dn = (LogNS(2,nSpecies) - LogNS(1,nSpecies))
!     if (dn >= 0) then
!        LogNS(0,nSpecies) = LogNS(1,nSpecies) - dn
!        LogNS(-1,nSpecies) = LogNS(0,nSpecies) - dn
!     else
!        LogNS(0,nSpecies) = LogNS(1,nSpecies) + dn
!        LogNS(-1,nSpecies) = LogNS(0,nSpecies) + dn
!     endif
!  endif

!  do iSpecies=nSpecies+1, nSpeciesTotal
!     dn = (LogNS(2,iSpecies) - LogNS(1,iSpecies))
!     LogNS(0,iSpecies) = LogNS(1,iSpecies) - dn
!     LogNS(-1,iSpecies) = LogNS(0,iSpecies) - dn
!!     if (LogNS( 0,iSpecies) < -30) LogNS( 0,iSpecies) = -30.0
!!     if (LogNS(-1,iSpecies) < -30) LogNS(-1,iSpecies) = -30.0
!  enddo

!  dn = (LogNS(2,iN_4S_) - LogNS(1,iN_4S_))
!  LogNS(0,iN_4S_) = LogNS(1,iN_4S_) - dn
!  LogNS(-1,iN_4S_) = LogNS(0,iN_4S_) - dn

  ! If you want to slip:
  !  Vel_GD(-1:0,iEast_)  = Vel_GD(1,iEast_)
  !  Vel_GD(-1:0,iNorth_) = Vel_GD(1,iNorth_)

  !-----------------------------------------------------------
  ! Top
  !-----------------------------------------------------------

  ! Slip flow at the top

  Vel_GD(nAlts+1:nAlts+2,iEast_)  = Vel_GD(nAlts,iEast_)
  Vel_GD(nAlts+1:nAlts+2,iNorth_) = Vel_GD(nAlts,iNorth_)

  IVel(nAlts+1:nAlts+2,iEast_)  = IVel(nAlts,iEast_)
  IVel(nAlts+1:nAlts+2,iNorth_) = IVel(nAlts,iNorth_)

  ! Things can go up or down in the ions

  IVel(nAlts+1,iUp_)   =  IVel(nAlts,iUp_)
  IVel(nAlts+2,iUp_)   =  IVel(nAlts-1,iUp_)

!  if(SZAVertical .GE. PI/2.) then
!     if(IVel(nAlts,iUp_) .GE. 0.0) then
!        IVel(nAlts,iUp_) = - IVel(nAlts,iUp_)
!        IVel(nAlts+1,iUp_) = - IVel(nAlts,iUp_)
!        IVel(nAlts+2,iUp_) = - IVel(nAlts,iUp_)
!        else
!           IVel(nAlts+1:nAlts+2,iUp_) = IVel(nAlts,iUp_)
!        endif
!     endif

  ! We only let stuff flow out in the neutrals

  if(Vel_GD(nAlts,iUp_)>0.)then
     Vel_GD(nAlts+1:nAlts+2,iUp_) = Vel_GD(nAlts,iUp_)
     VertVel(nAlts+1,:) = VertVel(nAlts,:)
     VertVel(nAlts+2,:) = VertVel(nAlts,:)
  else
     ! Vel_GD(nAlts+1:nAlts+2,iUp_) = 0.0 ! -Vel(nAlts)
     Vel_GD(nAlts+1,iUp_) = -Vel_GD(nAlts,iUp_)
     Vel_GD(nAlts+2,iUp_) = -Vel_GD(nAlts-1,iUp_)
     VertVel(nAlts+1,:) = -VertVel(nAlts,:)
     VertVel(nAlts+2,:) = -VertVel(nAlts-1,:)
  endif

  ! Constant temperature (zero gradient)

  Temp(nAlts+1) = Temp(nAlts)
  Temp(nAlts+2) = Temp(nAlts)

  dn = (LogRho(nAlts) - LogRho(nAlts-1))
  LogRho(nAlts+1) = LogRho(nAlts) + dn
  LogRho(nAlts+2) = LogRho(nAlts+1) + dn

  ! Limit the slope of the ion density

  Lst = mod(UTime/3600.0+Lon/15.0,24.0)

  do iSpecies=1,nIons-1 !Advect
     dn = (LogINS(nAlts,iSpecies) - LogINS(nAlts-1,iSpecies))

!     write(*,*) dn, dn * (2.0-cos(lst/12*Pi)), -0.01*LogINS(nAlts,iSpecies), lst
     fac = ((1.5-cos(lst/12*Pi))/1.75)**2.0
     dn = dn * fac
     if (dn > -0.01*LogINS(nAlts,iSpecies)) dn = -0.01*LogINS(nAlts,iSpecies)

!     dn = dn*(0.2*sin(Pi/2.-SZAVertical)+1.)
!     dn = 0.0

!     if (SZAVertical > Pi/2 .and. abs(lat) > 15.0 .and. abs(lat) < 70.0) then
!        dn = dn*abs(sin(lst*pi/12))*0.5
!        !write(*,*) lst,iSpecies, SZAVertical*180.0/3.14159, dn
!     endif

     LogINS(nAlts+1,iSpecies) = LogINS(nAlts,iSpecies) + dn
     LogINS(nAlts+2,iSpecies) = LogINS(nAlts+1,iSpecies) + dn
  enddo

  ! Hydrostatic pressure for the neutrals

  do iSpecies=1,nSpecies
     do iAlt = nAlts+1, nAlts+2
        InvScaleHeightS = -Gravity_G(iAlt) * &
             Mass(iSpecies) / (Temp(iAlt)*Boltzmanns_Constant)
        LogNS(iAlt,iSpecies) = &
             LogNS(iAlt-1,iSpecies) - dAlt_F(iAlt)*InvScaleHeightS
     enddo
     if (LogNS(nAlts+1,iSpecies) > 75.0 .or. &
          LogNS(nAlts+2,iSpecies) > 75.0) then
        write(*,*) "======> bcs : ", iSpecies, 1.0e-3/InvScaleHeightS, &
             Gravity_G(nAlts), Mass(iSpecies), Temp(nAlts), &
             LogNS(nAlts,iSpecies), LogNS(nAlts+1,iSpecies), &
             dAlt_F(nAlts), LogNS(nAlts+2,iSpecies)
     endif
  enddo

end subroutine set_vertical_bcs

