!\
! ------------------------------------------------------------
! set_boundary
! ------------------------------------------------------------
!/

subroutine set_vertical_bcs(LogRho,LogNS,Vel_GD,Temp, LogINS, iVel, VertVel)

  ! Fill in ghost cells at the top and bottom

  use ModSizeGitm, only: nAlts
  use ModPlanet, only: nSpecies, Mass, nIons, &
                       IsEarth, iN2_,iNO_, iN4S_, iO_
  use ModGITM, only: TempUnit, iEast_, iNorth_, iUp_, iProc
  use ModInputs
  use ModConstants
  use ModTime, only: UTime, iJulianDay,currenttime
  use ModVertical, only: Lat, Lon, Gravity_G, Altitude_G, dAlt_F, iLon1D, iLat1D,iBlock1d
  use ModEUV, only: sza
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
  real    :: HP
  integer :: ierror

  integer, dimension(25) :: sw

  !-----------------------------------------------------------
  ! Bottom
  !-----------------------------------------------------------

  ! Zero winds at the boundary
  !write(*,*) "HorizontalVelocityBC", HorizontalVelocityBC
  Vel_GD(-1:0,iEast_)  = HorizontalVelocityBC * cosd(lat)
  
  Vel_GD(-1:0,iNorth_) = 0.0

  Vel_GD(-1:0,iUp_)    = 0.0
  VertVel(-1:0,:)      = 0.0
  IVel(-1:0,iUp_)      = 0.0

  ! Constant gradient at bottom of atmosphere for ions:

  do iSpecies=1,nIons-1
     !dn = (LogINS(2,iSpecies) - LogINS(1,iSpecies))
     ! Set a zero gradient at the bottom
     dn = 0.0
     LogINS(0,iSpecies) = LogINS(1,iSpecies) - dn
     LogINS(-1,iSpecies) = LogINS(0,iSpecies) - dn
  enddo

  
  LogINS(0,nIons) = alog(sum(exp(LogINS(0,1:nIons-1))))
  LogINS(1,nIons) = alog(sum(exp(LogINS(1,1:nIons-1))))

  if (nSpecies >= iN4S_) then

     dn = (LogNS(2,iN4S_) - LogNS(1,iN4S_))
     if (dn >= 0) then
        LogNS(0,iN4S_) = LogNS(1,iN4S_) - dn
        LogNS(-1,iN4S_) = LogNS(0,iN4S_) - dn
     else
        LogNS(0,iN4S_) = LogNS(1,iN4S_) + dn
        LogNS(-1,iN4S_) = LogNS(0,iN4S_) + dn
     endif

     !Maybe just make this zero...
     if(VertVel(1,iN4S_) .gt. 0.0 )then
        ! Don't allow upwelling of N4S
        VertVel( 0,iN4S_) = 0.0
        VertVel(-1,iN4S_) = 0.0
     else
        VertVel( 0,iN4S_) =  0.0
        VertVel(-1,iN4S_) =  0.0
     endif
     
  endif
 
  if (nSpecies >= iO_) then

     !dn = (LogNS(2,iO_) - LogNS(1,iO_))
     !Make this a fixed bottom boundary
     !if (dn >= 0) then
     !   LogNS(0,iO_) = LogNS(1,iO_)
     !   LogNS(-1,iO_) = LogNS(0,iO_)
     !else
     !   LogNS(0,iO_) = LogNS(1,iO_)
     !   LogNS(-1,iO_) = LogNS(0,iO_)
     !endif

     !This may also need to be zero...
     if(VertVel(1,iO_) .gt. 0.0 )then
        ! Don't allow upwelling of O
        VertVel( 0,iO_) = 0.0
        VertVel(-1,iO_) = 0.0
     else
        VertVel( 0,iO_) =  0.0
        VertVel(-1,iO_) =  0.0
     endif
  
  endif

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
  IVel(nAlts+2,iUp_)   =  IVel(nAlts+1,iUp_)

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
  
  !BP: Add in the lower boundary condition for the temperature and/or a check
  !    for the temperature > 100.0  
  Temp(0) = 215
  Temp(-1) = 215

  dn = (LogRho(nAlts) - LogRho(nAlts-1))
  LogRho(nAlts+1) = LogRho(nAlts) + dn
  LogRho(nAlts+2) = LogRho(nAlts+1) + dn

  !if (iProc .eq. 1 .and. iLon1d .eq. 3 .and. iLat1d .eq. 4) then
    !write(*,*) "Pre-BCs (LogINS): ", LogINS(45:52,2)
    !write(*,*) "Pre-BCs (IDensityS): ", iDensityS(iLon1d,iLat1d,45:52,iBlock, 2)
    !write(*,*) "Pre-BCs (IVel): ", IVel(45:52,iUp_) 
  !endif

  ! Limit the slope of the ion density
  do iSpecies=1,nIons-1
     dn = (exp(LogINS(nAlts,iSpecies)) - exp(LogINS(nAlts-1,iSpecies)))
     
     if (dn > 0) dn = -0.01*exp(LogINS(nAlts,iSpecies))

     if (exp(LogINS(nAlts,iSpecies)) + dn < 0) then
       LogINS(nAlts+1,iSpecies) = -24
     else
       LogINS(nAlts+1,iSpecies) = alog(exp(LogINS(nAlts,iSpecies)) + dn)
     endif
     if (exp(LogINS(nAlts+1,iSpecies)) + dn < 0 .or. LogINS(nAlts+1,iSpecies) .eq. -24) then
       LogINS(nAlts+2,iSpecies) = -24
     else
       LogINS(nAlts+2,iSpecies) = alog(exp(LogINS(nAlts+1,iSpecies)) + dn)
     endif
  enddo
   
  !If using O+ nightside Venus ions in ghost cells...
  if (useNightSideIons .and. isVenus) then
    if (abs(sza(iLon1d, iLat1d, iBlock1d))*180.0/pi > 90.0) then
      LogINS(nAlts+1:nAlts+2,iOP_) = alog(5.0e9)
    endif 
  endif

  LogINS(nAlts+1,nIons) = alog(sum(exp(LogINS(nAlts+1,1:nIons-1))))
  LogINS(nAlts+2,nIons) = alog(sum(exp(LogINS(nAlts+2,1:nIons-1))))

  ! Hydrostatic pressure for the neutrals

  do iSpecies=1,nSpecies
     do iAlt = nAlts+1, nAlts+2
        InvScaleHeightS = -Gravity_G(iAlt) * &
             Mass(iSpecies) / (Temp(iAlt)*Boltzmanns_Constant)
        LogNS(iAlt,iSpecies) = &
             LogNS(iAlt-1,iSpecies) - dAlt_F(iAlt)*InvScaleHeightS
        if (LogNS(nAlts+1,iSpecies) > 75.0 .or. &
             LogNS(nAlts+2,iSpecies) > 75.0) then
           write(*,*) "======> bcs : ", iSpecies, 1.0e-3/InvScaleHeightS, &
                Gravity_G(nAlts), Mass(iSpecies), Temp(nAlts), &
                LogNS(nAlts,iSpecies), LogNS(nAlts+1,iSpecies), &
                dAlt_F(nAlts), LogNS(nAlts+2,iSpecies)
        endif
     enddo
  enddo

end subroutine set_vertical_bcs

