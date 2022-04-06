!--------------------------------------------------------------
!  Corrections by S. W. Bougher (9/28/07)
!  -- atomic oxygen replaces O2 in mean mass and scale height
!  -- artificially set all ion densities to 1.0 m-3 (not 1.0E+06 m-3)
!  Corrections by S. W. Bougher (11/07/12)
!  -- artificially set all ion densities to 1.0 m-3 (not 1.0E+04 m-3)
!  Corrections by S. W. Bougher (11/28/12)
!  -- artificially set O2+ and CO2+ ion densities to 1.0 m-3 (LBC at 80 km)
!--------------------------------------------------------------

subroutine get_msis_temperature(lon, lat, alt, t, h)

  use ModTime
  use ModInputs
  use ModPlanet
  use ModGITM

  implicit none

  real, intent(in) :: lon, lat, alt
  real, intent(out) :: t, h

! real :: nCO2, nO2, nN2, nCO, m, r, g
! real :: nCO2, nO,  nN2, nCO, m, r, g
  real :: nCO2, nOX,  nN2, nCO, nO2, m, r, g
  integer :: i

  i = 1
! do while (alt >= newalt(i))
  do while ( (alt >= newalt(i)) .and.  (i <= nAlts+2) )
     i = i + 1
  enddo
  i = i - 1

  t = InTemp(i)
  nCO2 = InNDensityS(i,iCO2_)
  nO2 = InNDensityS(i,iO2_)
! nO   = InNDensityS(i,iO_)
  nOX  = InNDensityS(i,iO_)
  nCO = InNDensityS(i,iCO_)
  nN2 = InNDensityS(i,iN2_)

! m = (nCO2 * mass(iCO2_) + &
!      nO2 * mass(iO2_) + &
!      nN2 * mass(iN2_) + &
!      nCO * mass(iCO_)) / (nCO2 + nO2 + nN2 + nCO)
! m = (nCO2 * mass(iCO2_) + &
!      nO   * mass(iO_) + &
!      nN2 * mass(iN2_) + &
!      nCO * mass(iCO_)) / (nCO2 + nO  + nN2 + nCO)

  m = (nCO2 * mass(iCO2_) + &
       nOX  * mass(iO_) + &
       nN2 * mass(iN2_) + &
       nCO * mass(iCO_) + &
       nO2 * mass(iO2_)) / (nCO2 + nOX + nN2 + nCO + nO2)

  r = RBody + alt
!  g = Gravitational_Constant * (RBody/r) ** 2
  g = Gravitational_Constant

  h = Boltzmanns_Constant * t / (m*g)

end subroutine get_msis_temperature


subroutine init_msis

  use ModPlanet
  use ModGITM
  use ModEUV
  use ModInputs

  implicit none

  integer, parameter:: ninitialAlts = 25
  integer :: iBlock
  integer :: iiLon,iiLat,iiAlt,iminiono,ialtlow(1), TimeArray(7),ilatlow(1)
  integer :: iLon,iLat,iAlt, iSpecies, iIon, iError, jlat,klon,iline

  logical :: Done = .False.,NotStarted = .True.
  character (len=iCharLen_) :: cLine
  real :: inDensities(9),altlow,althigh,latlow,lathigh
  real :: ralt, invAltDiff, altFind, altdiff, LogElectronDensity,dalt(nspeciestotal),alttemp(nInAlts)
  real, dimension(nInitialAlts) :: tempalt,LogInitialDensity,InitialEDensity,InitialAlt


  SurfaceAlbedo(:,:,:) = 0.0

  do iblock = 1, nblocks
     do ilon = 1,nlons
        do ilat = 1,nlats
           jlat = nint((latitude(ilat,iblock)*180.0/pi+93.75)/7.5)
           klon = nint((longitude(ilon,iblock)*180.0/pi+5.0)/10.0)
           SurfaceAlbedo(ilon,ilat,iblock) = dummyalbedo(jlat,klon)
           tinertia(ilon,ilat,iblock) = dummyti(jlat,klon)
      enddo
     enddo
  enddo

  if (DoRestart) return

  do iBlock = 1, nBlocks
     !write(*,*) "init_msis"
   do iLat = -1, nLats + 2
    iiLat = min(max(iLat,1),nLats)
    do iLon = -1, nLons + 2
     iiLon = min(max(iLon,1),nLons) 
     

     !BP (initialize all to 1)
     IDensityS(:,:,:,:,:) = 1.0e-24
     Temperature(:,:,:,:) = 150.0
  
     do iAlt = -1, nalts + 2
       alttemp = newalt

       altFind = altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0
       where(altfind - alttemp .lt. 0) alttemp = -1.0e9

       ialtlow = maxloc(alttemp)

       if (ialtlow(1) .eq. ninalts) ialtlow(1) = ialtlow(1) - 1

         altlow = newalt(ialtlow(1))
         althigh = newalt(ialtlow(1)+1)

         invaltdiff = 1/(althigh - altlow)

         if (altFind .lt. newalt(1)) then

           dalt = (altlow-altFind)*(InNDensityS(ialtlow(1) + 1,:) - &
                  inNDensityS(ialtlow(1),:)) * &
                  invAltDiff
           NDensityS(iLon,iLat,ialt,:,iBlock) = &
             inndensitys(ialtlow(1),:) - dalt
           ralt = (altlow-altFind)*(InTemp(ialtlow(1) + 1)-inTemp(ialtlow(1))) * &
                  invAltDiff
           Temperature(iLon,iLat,iAlt,iBlock) = InTemp(ialtlow(1)) - ralt

         else
           if (altFind .gt. newalt(ninalts)) then
             dalt = (altFind-althigh)* &
                    (inndensitys(ialtlow(1) + 1,:) - &
                    inndensitys(ialtlow(1),:)) * invAltDiff
             NDensityS(iLon,iLat,ialt,:,iBlock) = Inndensitys(ialtlow(1) + 1,:) + dalt
             ralt = (altFind-althigh)*(InTemp(ialtlow(1) + 1)-inTemp(ialtlow(1))) * &
                    invAltDiff
             Temperature(iLon,iLat,iAlt,iBlock)  = InTemp(ialtlow(1)) + ralt

           else 
             dalt = (Althigh - altFind)*(inNDensitys(ialtlow(1) + 1,:) - &
                    Inndensitys(ialtlow(1),:)) * invAltDiff
             NDensityS(iLon,iLat,ialt,:,iBlock) = inNDensityS(ialtlow(1) + 1,:) - dalt 
             ralt = (althigh-altFind)*(InTemp(ialtlow(1) + 1)-inTemp(ialtlow(1))) * &
                    invAltDiff
             Temperature(iLon,iLat,iAlt,iBlock) = InTemp(ialtlow(1)) - ralt
           endif
         endif
       enddo
      enddo! end iLon loop                                                                   
     enddo ! end iLat loop     


     !Initial Winds are -100 m/s (westward) unless otherwise specificed by the HorizontalBC
     do iLon = 1,nLons
       !write(*,*) "Latitude:", latitude(iLon,iBlock)
       Velocity(:,:,:,iEast_,iBlock)= HorizontalVelocityBC * cos(Latitude(iLon,iBlock))
     enddo
     !\
     ! Altitude Ghost Cells

     Temperature(:,:,-1,iBlock) = Temperature(:,:,1,iBlock)
     Temperature(:,:,0,iBlock) = Temperature(:,:,1,iBlock)

     Temperature(:,:,nAlts+1,iBlock) = Temperature(:,:,nAlts,iBlock)
     Temperature(:,:,nAlts+2,iBlock) = Temperature(:,:,nAlts,iBlock)
     eTemperature(:,:,nAlts+1,iBlock) = eTemperature(:,:,nAlts,iBlock)
     eTemperature(:,:,nAlts+2,iBlock) = eTemperature(:,:,nAlts,iBlock)

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

    !\
     ! Calculating MeanMajorMass -----------------------------
     !/

     !\
     ! Initialize MeanMajorMass to 0.0
     !/
     NDensityS = exp(nDensityS)

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


     TempUnit(-1:nLons+2,-1:nLats+2,-1:nAlts+2) = &
          MeanMajorMass(-1:nLons+2,-1:nLats+2,-1:nAlts+2)/&
          Boltzmanns_Constant

     !\
     ! Initialize Rho to 0.0
     !/

     Rho(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iBlock) = 0.0

     Temperature(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iBlock) = &
          Temperature(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iBlock) / &
          TempUnit(-1:nLons+2,-1:nLats+2,-1:nAlts+2)


     Rho(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iBlock) = &
          MeanMajorMass(-1:nLons+2,-1:nLats+2,-1:nAlts+2)* &
          NDensity(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iBlock)


     call calc_electron_temperature(iBlock)  

  enddo

end subroutine init_msis

subroutine msis_bcs(iJulianDay,UTime,Alt,Lat,Lon,Lst, &
             F107A,F107,AP,LogNS, Temp, LogRho)

  write(*,*) "You can not use MSIS with any planet except Earth!!!"
  write(*,*) "If you ARE running Earth, then make the code again, using"
  write(*,*) "configure Earth ; make"
  call stop_gitm("I can not continue...")

end subroutine msis_bcs


