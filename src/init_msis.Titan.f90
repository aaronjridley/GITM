!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine get_msis_temperature(lon, lat, alt, t, h)

  use ModTime
  use ModInputs
  use ModPlanet
  use ModGITM

  implicit none

  real, intent(in) :: lon, lat, alt
  real, intent(out) :: t, h

!! At Titan the only major species are CH4, N2, and H2 for the background (H2 is still pretty minor)

  real :: nN2, nCH4, nH2, m, r, g
  integer :: i

  i = 1
  do while ( (alt >= newalt(i)) .and.  (i <= nAlts+2) ) 
     i = i + 1
  enddo
  i = i - 1

  t = InTemp(i)
  nN2 = InNDensityS(i,iN2_)
  nCH4 = InNDensityS(i,iCH4_)
  nH2 = InNDensityS(i,iH2_)
  
  m = (nN2 * mass(iN2_) + &
       nCH4  * mass(iCH4_) + &
       nH2 * mass(iH2_) ) / (nN2 + nCH4 + nH2 )

  r = RBody + alt
  g = Gravitational_Constant * (RBody/r) ** 2
  h = Boltzmanns_Constant * t / (m*g)

end subroutine get_msis_temperature


subroutine init_msis

  use ModPlanet
  use ModGITM
  use ModEUV
  use ModInputs, only: iDebugLevel, iInputUnit_

  implicit none

  integer :: iBlock

  integer :: iiLon,iiLat,iiAlt
  integer :: iLon,iLat,iAlt, iSpecies, iIon

! Variables for the Asymmetric Wind Profile at the
! Lower Boundary

  real :: PeakLat
  real :: MinVelocity
  real :: RiseBeginLat, RiseEndLat
  real :: StandardT, RequiredT
  real :: RatioT
  real :: DelLat
  real :: NewLat
  real :: Freq

!! CH4 Mixing Ratio Fixes
  real :: DesiredCH4
  real :: TestCH4
  real :: CH4Scaling

!! CH4 Mixing Ratio Fixes
  real :: DesiredH2
  real :: TestH2
  real :: H2Scaling



  do iBlock = 1, nBlocks

     do iLat = -1, nLats + 2
        do iLon = -1, nLons + 2

           Temperature(iLon,iLat,-1:nAlts+2,iBlock) =  &
                InTemp(-1:nAlts+2)

           eTemperature(iLon,iLat,-1:nAlts+2,iBlock) =  &
                IneTemp(-1:nAlts+2)

           ITemperature(iLon,iLat,-1:nAlts+2,iBlock) =  &
                InTemp(-1:nAlts+2)

           NDensityS(iLon,iLat,-1:nAlts + 2,1:nSpeciesTotal,iBlock) =  1.0

           do iSpecies = 1, nSpeciesTotal

!! 0.42 Seemed too low
              NDensityS(iLon,iLat,-1:nAlts+2,iSpecies,iBlock) =  &
                   InNDensityS(-1:nAlts+2,iSpecies)*(1.0e+06)*(0.50)

           enddo ! end inner ispecies loop

           IDensityS(iLon,iLat,-1:nAlts + 2,ie_,iBlock) = 0.0

           do iIon = 1, nIons-1

             IDensityS(iLon,iLat,-1:nAlts + 2,iIon,iBlock) =  1.0e6

              IDensityS(iLon,iLat,-1:nAlts + 2,ie_,iBlock) = &
                   IDensityS(iLon,iLat,-1:nAlts + 2,ie_,iBlock) + &
                   IDensityS(iLon,iLat,-1:nAlts + 2,iIon,iBlock) 

           enddo ! end inner ispecies loop

        enddo! end iLon loop
     enddo ! end iLat loop

      PeakLat = 45.0*(pi/180.0)
      MinVelocity = 10.0

      RiseBeginLat = 30.0*(pi/180.0)
      RiseEndLat = 60.0*(pi/180.0)

      StandardT = 2.0*pi
      RequiredT = abs(RiseEndLat - RiseBeginLat)
      RatioT = StandardT/RequiredT

      DelLat = (pi/2.0 - RiseEndLat) 
      Freq = (pi/2.0)/DelLat

    do iLat = -1, nLats + 2
      do iLon = -1, nLons + 2


! Nothern Hemisphere

         if (Latitude(iLat,iBlock) .ge. 0.0) then
           if ( Latitude(iLat,iBlock) .le. RiseBeginLat) then
               Velocity(iLon,iLat,-1:1,iEast_,iBlock) = MinVelocity
           endif

           if ( (Latitude(iLat,iBlock) .ge. RiseBeginLat) .and. (Latitude(iLat,iBlock) .le. RiseEndLat) ) then
               Velocity(iLon,iLat,-1:1,iEast_,iBlock) = 30.0 - 20.0*cos(RatioT*Latitude(iLat,iBlock))
           endif

           if (Latitude(iLat,iBlock) .gt. RiseEndLat) then
               NewLat = Freq*(abs(Latitude(iLat,iBlock)) - RiseEndLat)
               Velocity(iLon,iLat,-1:1,iEast_,iBlock) = &
                      MinVelocity*( cos(NewLat) )**2.0
           endif
                
         endif  ! Northern Hemisphere

! Southern Hemisphere
         if (Latitude(iLat,iBlock) .lt. 0.0) then
               Velocity(iLon,iLat,-1:1,iEast_,iBlock) = MinVelocity
               if ( abs(Latitude(iLat,iBlock)) .gt. RiseEndLat) then

                    NewLat = Freq*(abs(Latitude(iLat,iBlock)) - RiseEndLat)

                    Velocity(iLon,iLat,-1:1,iEast_,iBlock) = &
                           MinVelocity*( cos(NewLat) )**2.0
               endif
         endif  ! Southern Hemisphere

      enddo !iLon
    enddo !iLat
      

!     where(NDensityS < 1.0e+03)
!        NDensityS = 1.0e+03
!     end where
!
!     where(IDensityS < 1.0e+03)
!        IDensityS = 1.0e+03
!     end where

     !\
     ! Altitude Ghost Cells

     Temperature(:,:,-1,iBlock) = Temperature(:,:,1,iBlock)
     Temperature(:,:,0,iBlock) = Temperature(:,:,1,iBlock)

     Temperature(:,:,nAlts+1,iBlock) = Temperature(:,:,nAlts,iBlock)
     Temperature(:,:,nAlts+2,iBlock) = Temperature(:,:,nAlts,iBlock)

     !\
     ! Longitude Ghost Cells
  
     Temperature(-1,:,:,iBlock) = Temperature(1,:,:,iBlock)
     Temperature(0,:,:,iBlock) = Temperature(1,:,:,iBlock)

     Temperature(nLons+1,:,:,iBlock) = Temperature(nLons,:,:,iBlock)
     Temperature(nLons+2,:,:,iBlock) = Temperature(nLons,:,:,iBlock)

     !\
     ! Latitude Ghost Cells
     
     Temperature(:,-1,:,iBlock) = Temperature(:,1,:,iBlock)
     Temperature(:,0,:,iBlock) = Temperature(:,1,:,iBlock)

     Temperature(:,nLats+1,:,iBlock) = Temperature(:,nLats,:,iBlock)
     Temperature(:,nLats+2,:,iBlock) = Temperature(:,nLats,:,iBlock)


     MeanMajorMass(-1:nLons+2,-1:nLats+2,-1:nAlts+2) = 0.0
     MeanIonMass(-1:nLons+2,-1:nLats+2,-1:nAlts+2) = 0.0


! Calculate MeanMajorMass -----------------------------
! Calculate TempUnit -----------------------------

     do iLat = -1,nLats + 2
        do iLon = -1,nLons + 2
           do iAlt = -1,nAlts + 2

              NDensity(iLon,iLat,iAlt,iBlock) = 0.0
           
              do iSpecies = 1,nSpeciesTotal
                 NDensity(iLon,iLat,iAlt,iBlock) = &
                      NDensity(iLon,iLat,iAlt,iBlock) + &
                      NDensityS(iLon,iLat,iAlt,iSpecies,iBlock)
              enddo

              do iSpecies = 1,nSpeciesTotal
                 MeanMajorMass(iLon,iLat,iAlt) = &
                      MeanMajorMass(iLon,iLat,iAlt) + &
                      Mass(iSpecies)*NDensityS(iLon,iLat,iAlt,iSpecies,iBlock)/ &
                      NDensity(iLon,iLat,iAlt,iBlock) 
              enddo

              do iIon = 1,nIons - 1
                 MeanIonMass(iLon,iLat,iAlt) = &
                      MeanIonMass(iLon,iLat,iAlt) + &
                      MassI(iIon)*IDensityS(iLon,iLat,iAlt,iIon,iBlock)/ &
                      IDensityS(iLon,iLat,iAlt,ie_,iBlock) 
              enddo


           enddo
        enddo
     enddo

     TestCH4 = NDensityS(1,1,1,iCH4_,iBlock)/NDensity(1,1,1,iBlock)
     DesiredCH4 = 0.014
     CH4Scaling = DesiredCH4/TestCH4

     do iLat = -1,nLats + 2
        do iLon = -1,nLons + 2
           do iAlt = -1,nAlts + 2

              NDensityS(iLon,iLat,iAlt,iCH4_,iBlock) = &
                  CH4Scaling*NDensityS(iLon,iLat,iAlt,iCH4_,iBlock)

              NDensityS(iLon,iLat,iAlt,i13CH4_,iBlock) = &
              (1.0/79.0)*NDensityS(iLon,iLat,iAlt,iCH4_,iBlock)

              NDensityS(iLon,iLat,iAlt,i15N2_,iBlock) = &
              (2.0/140.0)*NDensityS(iLon,iLat,iAlt,iN2_,iBlock)

           enddo
        enddo
     enddo

     TestH2 = NDensityS(1,1,1,iH2_,iBlock)/NDensity(1,1,1,iBlock)
     DesiredH2 = 0.00400
     H2Scaling = DesiredH2/TestH2

     do iLat = -1,nLats + 2
        do iLon = -1,nLons + 2
           do iAlt = -1,nAlts + 2

              NDensityS(iLon,iLat,iAlt,iH2_,iBlock) = &
                  H2Scaling*NDensityS(iLon,iLat,iAlt,iH2_,iBlock)

           enddo
        enddo
     enddo

     do iLat = -1,nLats + 2
        do iLon = -1,nLons + 2

              NDensityS(iLon,iLat,-1:5,iAr_,iBlock) = &
                  (4.32e-05)*NDensity(iLon,iLat,-1:5,iBlock)

              NDensityS(iLon,iLat,6:nAlts,iAr_,iBlock) = &
                  (1.00e-10)*NDensity(iLon,iLat,6:nAlts,iBlock)

        enddo
     enddo

     TempUnit(-1:nLons+2,-1:nLats+2,-1:nAlts+2) = &
          MeanMajorMass(-1:nLons+2,-1:nLats+2,-1:nAlts+2)/&
          Boltzmanns_Constant

     Rho(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iBlock) = 0.0

     Temperature(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iBlock) = &
          Temperature(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iBlock) / &
          TempUnit(-1:nLons+2,-1:nLats+2,-1:nAlts+2)

     Rho(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iBlock) = &
          MeanMajorMass(-1:nLons+2,-1:nLats+2,-1:nAlts+2)* &
          NDensity(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iBlock)

  enddo

end subroutine init_msis

subroutine msis_bcs(iJulianDay,UTime,Alt,Lat,Lon,Lst, &
             F107A,F107,AP,LogNS, Temp, LogRho)

  write(*,*) "You can not use MSIS with any planet except Earth!!!"
  write(*,*) "If you ARE running Earth, then make the code again, using"
  write(*,*) "configure Earth ; make"
  call stop_gitm("I can not continue...")

end subroutine msis_bcs
