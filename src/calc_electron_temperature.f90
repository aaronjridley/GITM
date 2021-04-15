
subroutine calc_electron_temperature(iBlock)

  !  Take Tion and Telec values from empirical datasets from Viking (settei.F)
  !  -- Pre-MAVEN
  !  Take Telec values from Ergun et al (2015) plus Linear Fit to Tn at 130 km
  !  -- Post-MAVEN  (dayside low SZA experiment for dayside orbit application)

  use ModGITM
  use ModPlanet, only : ialtminiono

  implicit none

  integer, intent(in) :: iBlock
  integer :: iLon,iLat,iAlt,iminiono,k130
  real :: Alt,TN130
! real,parameter :: TL = 510.
  real,parameter :: TL2 = 510.
  real,parameter :: TH = 3140.
  real,parameter :: Z0 = 241.

  do iLon = 1, nLons
    do iLat = 1, nLats
      do iAlt = 1,nAlts
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







