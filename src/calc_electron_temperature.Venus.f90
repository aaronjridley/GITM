
subroutine calc_electron_temperature(iBlock)
  use ModGITM
  use ModPlanet, only : ialtminiono

  implicit none

  integer, intent(in) :: iBlock
  integer :: iLon,iLat,iAlt

  do iLon = -1, nLons+2
    do iLat = -1, nLats+2
      do iAlt = -1,nAlts+2
        iTemperature(iLon,iLat,iAlt,iBlock) = Temperature(iLon,iLat,iAlt,iBlock)*TempUnit(iLon,iLat,iAlt)
        if (Altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0 > 140.0) then
          eTemperature(iLon,iLat,iAlt,iBlock) = 10**(3.471 - 1921.9/(Altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0 &
                                                - 98.078)**2. + &
                                                8.5257*Altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0*1.0E-04)
        else
          eTemperature(iLon,iLat,iAlt,iBlock) = Temperature(iLon,iLat,iAlt,iBlock)*TempUnit(iLon,iLat,iAlt)
        endif
      enddo
    enddo
  enddo
end subroutine calc_electron_temperature







