!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine calc_electron_temperature(iBlock)

  !  Take values from empirical datasets fropm Viking (settei.F)

  use ModGITM
  use ModPlanet, only : ialtminiono

  implicit none

  integer, intent(in) :: iBlock
  integer :: iLon,iLat,iAlt,iminiono
  real :: Alt

  call report("Electron Density", 2)
  !Electron Temperature
  !FOX[93] FORMULATION FOR DY TE FROM ROHRBAUGH ET. AL. [79]

  do iLon = 1, nLons
     do iLat = 1, nLats
        iMinIono = iAltMinIono(iLon,iLat,iBlock)
        do ialt = iminiono, nAlts

           Alt = Altitude_GB(iLon,iLat,iAlt,iBlock) /1000.0
           if (Alt < 130.0) eTemperature(iLon,iLat,iAlt,iBlock) = &
                Temperature(iLon,iLat,iAlt,iBlock) * TempUnit(iLon,iLat,iAlt)

           if (Alt > 180.0) eTemperature(iLon,iLat,iAlt,iBlock) = &
                4200.0 - 3750.0*exp((180-Alt)/89.6)

           if (Alt >= 130.0 .and. Alt <= 180.0) eTemperature(iLon,iLat,iAlt,iBlock) = &
                700.0-536.0*exp((130.0 - Alt)/65.4)

        enddo
     enddo
  enddo

  !Ion Temperature
  !FOX[93] FORMULATION FOR DY TE FROM ROHRBAUGH ET. AL. [79]

  do iLon = 1, nLons
     do iLat = 1, nLats
        iMinIono = iAltMinIono(iLon,iLat,iBlock)
        do ialt = iminiono, nAlts

           Alt = Altitude_GB(iLon,iLat,iAlt,iBlock) /1000.0
           if (Alt < 180.0) iTemperature(iLon,iLat,iAlt,iBlock) = &
                Temperature(iLon,iLat,iAlt,iBlock) * TempUnit(iLon,iLat,iAlt)

           if (Alt > 300.0) iTemperature(iLon,iLat,iAlt,iBlock) = &
                eTemperature(iLon,iLat,iAlt,iBlock)

           if (Alt >= 180.0 .and. Alt <= 300.0) iTemperature(iLon,iLat,iAlt,iBlock) = &
                10.0**(2.243+(Alt - 180.0)/95.0)

        enddo
     enddo
  enddo

end subroutine calc_electron_temperature







