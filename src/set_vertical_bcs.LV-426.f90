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
  use ModPlanet, only: nSpecies, nSpeciesTotal, nIonsAdvect, Mass, nIons
  use ModGITM, only: iEast_, iNorth_, iUp_
  use ModInputs
  use ModConstants
  use ModVertical, only: &
       Lat, Lon, Gravity_G, Altitude_G, dAlt_F, iLon1D, iLat1D, iBlock1D

  implicit none

  real, intent(inout) :: &
       LogRho(-1:nAlts+2), &
       LogNS(-1:nAlts+2,nSpecies), &
       LogINS(-1:nAlts+2,nIonsAdvect), &
       Vel_GD(-1:nAlts+2,3), &
       IVel(-1:nAlts+2,3), &
       Temp(-1:nAlts+2), &
       VertVel(-1:nAlts+2,nSpecies)

  integer :: iSpecies, iAlt
  real    :: InvScaleHeightS, InvScaleHgt, dn
  integer :: ierror

  !-----------------------------------------------------------
  ! Bottom
  !-----------------------------------------------------------

  ! Don't Let the winds blow
  Vel_GD(-1:0,iEast_)  = 0.0
  Vel_GD(-1:0,iNorth_) = 0.0
  Vel_GD(0,iUp_)       = -Vel_GD(0,iUp_)
  Vel_GD(-1,iUp_)      = -Vel_GD(1,iUp_)
  VertVel(0,:)         = -VertVel(0,:)
  VertVel(-1,:)        = -VertVel(1,:)
  IVel(-1:0,iUp_)      = 0.0

  Temp(0)       = Temp(1)
  Temp(-1)      = Temp(0)

  do iSpecies=1, nSpeciesTotal
     dn = (LogNS(2,iSpecies) - LogNS(1,iSpecies))
     LogNS(0,iSpecies) = LogNS(1,iSpecies) - dn
     LogNS(-1,iSpecies) = LogNS(0,iSpecies) - dn
  enddo

  do iSpecies=1,nIonsAdvect
     dn = (LogINS(2,iSpecies) - LogINS(1,iSpecies))
     LogINS(0,iSpecies) = LogINS(1,iSpecies) - dn
     LogINS(-1,iSpecies) = LogINS(0,iSpecies) - dn
  enddo

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

  do iSpecies=1,nIonsAdvect
     dn = (LogINS(nAlts,iSpecies) - LogINS(nAlts-1,iSpecies))
!     if (dn < 0.75*LogINS(nAlts,iSpecies) .and. dn > 0) &
!          dn = 0.75*LogINS(nAlts,iSpecies)
     if (dn > -0.25*LogINS(nAlts,iSpecies)) dn = -0.25*LogINS(nAlts,iSpecies)
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

