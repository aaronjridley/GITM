!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine calc_physics(iBlock)

  use ModGITM
  use ModTime
  use ModConstants
  use ModEUV
  use ModPlanet
  use ModInputs
  implicit none

  integer, external :: jday

  integer, intent(in) :: iBlock

  integer :: itime(7)
  real*8 :: VernalDay
  real :: JulianDayEq, DaysInYear, OrbitAngle
  real :: SinDec, CosDec

  real :: DayNumber,TimeHour
  real :: MeanAnomaly,DMeanAnomaly, TrueAnomaly

  real :: DNoverCY, tol, eccentricityDeg, dm
  real :: semimajoraxis, eccentricity,inclination,meanLongitude,longitudePerihelion
  real :: longitudeNode, argPerihelion, EccAnomaly, dEccAnomaly
  real :: heliocentricX, heliocentricY, heliocentricZ
  real :: x_ecl, y_ecl

  integer :: i
  integer :: iLon, iLat, iAlt

  call report("calc_physics",2)

  !\
  ! Solar Things
  !/

  !!!! Calculate orbital parameters based on E Standish, Solar System Dynamics, JPL,.
  !!!! No constant orbital speed assumption
  TimeHour = iTimeArray(4)+iTimeArray(5)/60.0+iTimeArray(6)/3600.0

  DayNumber = 367*iTimeArray(1)-7*(iTimeArray(1)+(iTimeArray(2)+9)/12)/4 &
     +275*iTimeArray(2)/9+iTimeArray(3)-730531.5+TimeHour/24.0

  DNoverCY = DayNumber/36525.0

  !Compute Keplerian elements at current date

  semimajoraxis       = semimajoraxis_0+semimajordot*DNoverCY     !AU
  eccentricity        = eccentricity_0+eccentricitydot*DNoverCY   !Rad
  inclination         = inclination_0+inclinationdot*DNoverCY     !deg
  meanLongitude       = meanLongitude_0+meanLongitudedot*DNoverCY !deg
  longitudeNode       = semimajoraxis_0+semimajordot*DNoverCY     !deg
  longitudePerihelion = longitudePerihelion_0+longitudePeriheliondot*DNoverCY !deg

  !Compute argument of perihelion and mean anomaly
  argPerihelion = longitudePerihelion - longitudeNode
  !computation of M for Jupiter and out is supposed to be modified by additional
  !terms for the time interval 3000BC - 30000AD... This probably doesn't matter.
  MeanAnomaly = meanLongitude - longitudePerihelion
  MeanAnomaly = mod(MeanAnomaly,360.0)
  if (MeanAnomaly > 180) MeanAnomaly = MeanAnomaly - 360

  !Need to solve Kepler's equation by iterating
  DEccAnomaly = 10000.0
  tol = 1.0e-6
  eccentricityDeg = (180/pi)*eccentricity
  EccAnomaly = MeanAnomaly+eccentricityDeg*sin(MeanAnomaly*pi/180)
  i = 0
  do while (abs(DEccAnomaly) > tol .and. i < 100)
    dM = MeanAnomaly - (eccAnomaly - eccentricityDeg*sin(pi/180*eccAnomaly))
    DEccAnomaly = dM/(1-eccentricity*cos(pi/180*eccAnomaly))
    eccAnomaly = eccAnomaly+dEccAnomaly
    i = i+1
  end do

  !Get heliocentric coordinates, TrueAnomaly and sunorbiteccentricity
  heliocentricX = semimajoraxis*(cos(EccAnomaly*pi/180.)-eccentricity)
  heliocentricY = semimajoraxis*(1-eccentricity**2)**(0.5)*sin(eccAnomaly*pi/180.)
  heliocentricZ = 0

  TrueAnomaly = atan2(heliocentricY,heliocentricX)*180/pi
  sunorbiteccentricity = sqrt(heliocentricX**2+heliocentricY**2)

  !convert to J2000 coordinates with x-axis aligned with vernal equinox so
  !we can get solar longitude in the coorect system.  We don't need z.
  x_ecl = heliocentricX * &
    (cos(argPerihelion*pi/180.)*cos(longitudeNode*pi/180.) - &
    sin(argPerihelion*pi/180.)*sin(longitudeNode*pi/180.)*cos(inclination*pi/180.)) + &
    heliocentricY * &
    (-sin(argPerihelion*pi/180.)*cos(longitudeNode*pi/180.) - &
    cos(argPerihelion*pi/180.)*sin(longitudeNode*pi/180.)*cos(inclination*pi/180.))

  y_ecl = heliocentricX * &
    (cos(argPerihelion*pi/180.)*sin(longitudeNode*pi/180.) + &
    sin(argPerihelion*pi/180.)*cos(longitudeNode*pi/180.)*cos(inclination*pi/180.)) + &
    heliocentricY * &
    (-sin(argPerihelion*pi/180.)*sin(longitudeNode*pi/180.) + &
    cos(argPerihelion*pi/180.)*cos(longitudeNode*pi/180.)*cos(inclination*pi/180.))

  !Calculate orbit angle, aka Ls. In this CS, need the angle from -x axis.
  orbitAngle = atan(y_ecl/x_ecl)
  if (x_ecl > 0) orbitangle = orbitangle+pi
  if (x_ecl < 0 .and. y_ecl > 0) orbitangle = orbitangle + 2*pi
  SunDeclination = atan(tan(Tilt*pi/180.)*sin(OrbitAngle))

  SinDec = sin(SunDeclination)
  CosDec = cos(SunDeclination)

  ! Updated to work at all Planets
  ! We need the equivalent number of "hours" per planet rotation,
  ! assuming that there are 24.0 LT hours per planet day
  LocalTime = mod((UTime/(Rotation_Period/24.0) + &
       Longitude(:,iBlock) * 24.0/TwoPi), 24.0)

  if (UseApex) &
       call SUBSOLR(iTimeArray(1),iJulianDay,iTimeArray(4),&
       iTimeArray(5),iTimeArray(6),SubsolarLatitude, &
       SubsolarLongitude)

  do iAlt=-1,nAlts+2

     !
     ! Compute Magnetic Local Time
     !

     if (UseApex .and. IsEarth) then
        do iLat=-1,nLats+2
           do iLon=-1,nLons+2
              call magloctm( &
                   MLongitude(iLon,iLat,iAlt,iBlock), &
                   SubsolarLatitude,   &
                   SubsolarLongitude,  &
                   MagneticPoleColat, &
                   MagneticPoleLon,   &
                   MLT(iLon,iLat,iAlt))
              if (mlt(iLon,iLat,iAlt) < 0) &
                   mlt(iLon,iLat,iAlt) = mlt(iLon,iLat,iAlt) + 24.0
           enddo
        enddo
     else

        do iLat=-1,nLats+2
           MLT(:,iLat,iAlt) = LocalTime
        enddo

        ! Since we go over the pole,
        !we have to move the points to the proper location:

        if (Latitude(0,iBlock) < -pi/2.0) then
           MLT(:,-1,iAlt)      = MLT(:,2,iAlt) + 12.0
           MLT(:,0,iAlt)       = MLT(:,1,iAlt) + 12.0
        endif

        if (Latitude(0,iBlock) > pi/2.0) then
           MLT(:,nLats+1,iAlt) = MLT(:,nLats,iAlt) + 12.0
           MLT(:,nLats+2,iAlt) = MLT(:,nLats-1,iAlt) + 12.0
        endif

     endif

  enddo

!  write(*,*) mlt(1,1,1)

  where (MLT > 12.0)
     MLT = MLT - 24.0
  end where

  do iLon = 1, nLons
     do iLat = 1, nLats
        sza(iLon, iLat,iBlock) =  &
             acos(SinDec*sin(Latitude(iLat,iBlock)) + &
             CosDec*CosLatitude(iLat,iBlock) * &
             cos(pi*(LocalTime(iLon)- 24.0/2.0)/(24.0/2.)))

        if (DtLTERadiation < Rotation_Period) then
           call calc_avesza(iLon,iLat,iBlock, SinDec, CosDec)
        endif

     enddo
  enddo

  call calc_scaled_euv

  do iAlt = 1, nAlts
     xSolar(:,:,iAlt) = RadialDistance_GB(1:nLons,1:nLats,iAlt,iBlock) &
          * cos(SZA(:,:,iBlock))
     ySolar(:,:,iAlt) = RadialDistance_GB(1:nLons,1:nLats,iAlt,iBlock) &
          * sin(SZA(:,:,iBlock))
  enddo

end subroutine calc_physics

!-----------------------------------------------------------------------------
! get_subsolar: A routine to get the latitude and longitude of the subsolar
!               point for a specified time and Vernal Equinox.
!
! Author: Alexey Morozov, UMich, December 2012
!
! Comments: Tested by Angeline G Burrell on December 26, 2012
!----------------------------------------------------------------------------

subroutine get_subsolar(CurrentTime, VernalTime, lon_sp, lat_sp)

  use ModConstants, only : pi
  use ModPlanet, only : Tilt, SecondsPerYear, Rotation_Period

  implicit none

  real*8, intent(in) :: CurrentTime, VernalTime
  real*8, intent(out) :: lon_sp, lat_sp

  lon_sp=(pi - ((CurrentTime - VernalTime)/Rotation_Period &
       - floor((CurrentTime - VernalTime)/Rotation_Period))*2*pi)
  if (lon_sp<0.) lon_sp=lon_sp+2*pi

  lat_sp=atan(tan(Tilt*pi/180.)*sin(2.*pi*(CurrentTime - VernalTime) &
       /SecondsPerYear))

end subroutine get_subsolar

!-----------------------------------------------------------------------------
! get_sza: A routine to get the sza given latitude and longitude
!
! Author: Ankit Goel, UMich, 21 April 2015
!
! Comments:
!----------------------------------------------------------------------------

subroutine get_sza(lon, lat, sza_calc)
  use ModGITM
  use ModTime
  use ModConstants
  use ModEUV
  use ModPlanet
  use ModInputs
  implicit none

  double precision, intent(in) :: lon, lat !lat, lon need to be in radians
  double precision, intent(out) :: sza_calc

  real(8) :: OrbitAngle_1, LocalTime_1
  real(8) :: SunDeclination_1, SinDec, CosDec

  OrbitAngle_1 = 2.*pi*(CurrentTime - VernalTime)/SecondsPerYear
  SunDeclination_1 = atan(tan(Tilt*pi/180.)*sin(OrbitAngle_1))

  SinDec = sin(SunDeclination_1)
  CosDec = cos(SunDeclination_1)

  LocalTime_1 = mod((UTime/3600.0 + &
       lon * HoursPerDay / TwoPi), HoursPerDay)

  sza_calc = acos(SinDec*sin(Lat) + &
                  CosDec*Cos(lat) * &
                  cos(pi*(LocalTime_1-HoursPerDay/2)/(HoursPerDay/2)))

end subroutine get_sza
