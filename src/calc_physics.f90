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
  real :: SunDeclination, SinDec, CosDec

  real :: DayNumber,TimeHour,Cy
  real :: MeanAnomaly,RMeanAnomaly,DMeanAnomaly
  real :: EccAnomaly1,EccAnomaly2,REccAnomaly1,REccAnomlay2
  real :: TrueAnomaly,RTrueAnomaly
  
  integer :: i
  integer :: iLon, iLat, iAlt

  call report("calc_physics",2)

  !\
  ! Solar Things
  !/

!\
! Orbital Elements procedure taken from 
! http://home.att.net/~srschmitt/planetorbits.html
!/

  TimeHour = iTimeArray(4)+iTimeArray(5)/60.0+iTimeArray(6)/3600.0
  
  DayNumber = 367*iTimeArray(1)-7*(iTimeArray(1)+(iTimeArray(2)+9)/12)/4 &
	   +275*iTimeArray(2)/9+iTimeArray(3)-730531.5+TimeHour/24.0
 
  CY = DayNumber/36525.0
  
  MeanAnomaly =(SunOrbit_D+((SunOrbit_E*CY)/3600.0)) - SunOrbit_C
  RMeanAnomaly = pi/180.0*MeanAnomaly
  DMeanAnomaly = RMeanAnomaly/(2.0*pi)
  
  if (DMeanAnomaly > 0.0) then
	 RMeanAnomaly = (2.0*pi)*(DMeanAnomaly-floor(DMeanAnomaly))
  else
	 RMeanAnomaly = (2.0*pi)*(DMeanAnomaly-ceiling(DMeanAnomaly))
  endif
  
  if (RMeanAnomaly < 0.0) then
	 RMeanAnomaly = (2.0*pi)+RMeanAnomaly
  endif
  
  EccAnomaly1 = RMeanAnomaly + SunOrbit_B*sin(RMeanAnomaly)*(1.0+ &
	 SunOrbit_B*cos(RMeanAnomaly))
  
  !  Calculate true anomaly  
  do i=1,10
	 EccAnomaly2 = EccAnomaly1
	 
	 EccAnomaly1 = EccAnomaly2 - ((EccAnomaly2-SunOrbit_B* &
		sin(EccAnomaly2)-RMeanAnomaly)/(1.0-SunOrbit_B* &
		cos(EccAnomaly2)))

  enddo

  TrueAnomaly = 2.0*atan(sqrt((1.0+SunOrbit_B)/ &
	 (1.0-SunOrbit_B))*tan(0.5*EccAnomaly1))
 
  if (TrueAnomaly < 0) then
	 TrueAnomaly = TrueAnomaly + (2.0*pi)
  endif
  
 !!!!!!! General SunOrbitEccentricity 
  SunOrbitEccentricity = SunOrbit_A*(1.0-SunOrbit_B**2.0)/(1.0+ &
	   SunOrbit_B*cos(TrueAnomaly))

 OrbitAngle = 2.*pi*(CurrentTime - VernalTime)/SecondsPerYear

 !!!!!!! Old SunOrbitEccentricity
!  SunOrbitEccentricity = SunOrbit_A                     + &
!                         SunOrbit_B*cos(OrbitAngle)    + &
!                         SunOrbit_C*sin(OrbitAngle)    + &
!                         SunOrbit_D*cos(2.*OrbitAngle) + &
!                         SunOrbit_E*sin(2.*OrbitAngle)


  SunDeclination = atan(tan(Tilt*pi/180.)*sin(OrbitAngle))

  SinDec = sin(SunDeclination)
  CosDec = cos(SunDeclination)

  LocalTime = mod((UTime/3600.0 + &
       Longitude(:,iBlock) * HoursPerDay / TwoPi), HoursPerDay)

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
             cos(pi*(LocalTime(iLon)-HoursPerDay/2)/(HoursPerDay/2)))

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
