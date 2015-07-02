!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine get_msis_temperature(lon, lat, alt, t, h)

  use ModIndicesInterfaces
  use ModTime
  use ModInputs
  use ModPlanet
  use ModGITM

  use EUA_ModMsis90, only: meter6, gtd6

  implicit none

  real, intent(in) :: lon, lat, alt
  real, intent(out) :: t, h

  real :: m, r, g

  !-------------------------------------------------------
  
  t = 100.0

  m = Mass(iH2_)

  r = RBody + alt
  g = Gravitational_Constant * (RBody/r) ** 2

  h = Boltzmanns_Constant * t / (m*g)

end subroutine get_msis_temperature

!--------------------------------------------------------------
!
!--------------------------------------------------------------

subroutine init_msis

  use ModGITM
  use ModInputs
  use ModConstants
  use ModPlanet
  use ModTime

  implicit none

  integer :: iBlock, iAlt, iLat, iLon, iSpecies, iyd
  real :: geo_lat, geo_lon, geo_alt

  if (DoRestart) return

  ! Initialize data

  do iBlock = 1, nBlocks
     do iAlt = -1, nAlts+2
        do iLon=-1,nLons+2
           do iLat=-1,nLats+2

              NDensityS(iLon,iLat,iAlt,iH2_,iBlock) = 0.0

              MeanMajorMass(iLon,iLat,iAlt) = Mass(iH2_)
  
              TempUnit(iLon,iLat,iAlt) = &
                   MeanMajorMass(iLon,iLat,iAlt)/ Boltzmanns_Constant

              Temperature(iLon,iLat,iAlt,iBlock) = 100.0

              Rho(iLon,iLat,iAlt,iBlock) = &
                   MeanMajorMass(iLon,iLat,iAlt) * &
                   NDensityS(iLon,iLat,iAlt,iH2_,iBlock)

              LogNS(iLon,iLat,iAlt,:,iBlock) = &
                   log(NDensityS(iLon,iLat,iAlt,:,iBlock))

              NDensity(iLon,iLat,iAlt,iBlock) = &
                   sum(NDensityS(iLon,iLat,iAlt,1:nSpecies,iBlock))

              Velocity(iLon,iLat,iAlt,:,iBlock) = 0.0

           enddo
        enddo
     enddo

     Rho(:,:,:,iBlock) = 0.0
     NDensity(:,:,:,iBlock) = 0.0

     do iSpecies=1,nSpecies

        NDensity(:,:,:,iBlock) = NDensity(:,:,:,iBlock) + &
             NDensityS(:,:,:,iSpecies,iBlock)

        Rho(:,:,:,iBlock) = Rho(:,:,:,iBlock) + &
             Mass(iSpecies)*NDensityS(:,:,:,iSpecies,iBlock)
        
     enddo

  enddo
 
end subroutine init_msis

!--------------------------------------------------------------
!
!--------------------------------------------------------------

subroutine msis_bcs(iJulianDay,UTime,Alt,Lat,Lon,Lst, &
     F107A,F107,AP,LogNS, Temp, LogRho, v)

  use ModTime, only : iTimeArray
  use ModPlanet

  implicit none

  integer, intent(in) :: iJulianDay
  real, intent(in) :: uTime, Alt, Lat, Lon, LST, f107a, f107
  real, intent(in):: ap
  real, intent(out) :: LogNS(nSpecies), Temp, LogRho, v(2)
  real :: h2

  h2 = 1.0e10
  LogNS(iH2_)  = alog(h2)

  Temp        = 100.0
  LogRho      = alog(h2*mass(iH2_))

  V(1) = 0.0
  V(2) = 0.0

end subroutine msis_bcs

