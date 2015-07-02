!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!--------------------------------------------------------------
!
!--------------------------------------------------------------

subroutine get_msis_temperature(lon, lat, alt, t, h)

  use ModTime
  use ModInputs
  use ModPlanet
  use ModGITM

  implicit none

  real, intent(in) :: lon, lat, alt
  real, intent(out) :: t, h

  real :: nCO2, nO2, nN2, nCO, m, r, g
  integer :: i

  i = 1
  do while (alt >= newalt(i)) 
     i = i + 1
  enddo
  i = i - 1

  t = InTemp(i)
  nCO2 = InNDensityS(i,iCO2_)
  nO2 = InNDensityS(i,iO2_)
  nCO = InNDensityS(i,iCO_)
  nN2 = InNDensityS(i,iN2_)
  
  m = (nCO2 * mass(iCO2_) + &
       nO2 * mass(iO2_) + &
       nN2 * mass(iN2_) + &
       nCO * mass(iCO_)) / (nCO2 + nO2 + nN2 + nCO)

  r = RBody + alt
!  g = Gravitational_Constant * (RBody/r) ** 2
  g = Gravitational_Constant

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

  write(*,*) "here!"

  do iBlock = 1, nBlocks

     write(*,*) '==> Now Initializing Mars Background Composition', iBlock   

!\
! Initializes the Planet with the same Chemistry as Above 
!/

     do iLat = -1, nLats + 2
        do iLon = -1, nLons + 2

           Temperature(iLon,iLat,-1:nAlts+2,iBlock) =  &
                InTemp(-1:nAlts+2)

           eTemperature(iLon,iLat,-1:nAlts+2,iBlock) =  &
                IneTemp(-1:nAlts+2)

           ITemperature(iLon,iLat,-1:nAlts+2,iBlock) =  &
                InITemp(-1:nAlts+2)

           NDensityS(iLon,iLat,-1:nAlts + 2,1:nSpeciesTotal,iBlock) =  1.0

           do iSpecies = 1, nSpeciesTotal

              ! Converts from cm^-3 to m^-3
              NDensityS(iLon,iLat,-1:nAlts+2,iSpecies,iBlock) =  &
                   InNDensityS(-1:nAlts+2,iSpecies)*(1.0e+06)

           enddo ! end inner ispecies loop

           !\
           ! These first few are from the input file read in earlier
           !

        
           !\
           ! This just arbitrarily sets the ionospheric densities to 1.0
           !

           do iIon = 1, nIons-1

              IDensityS(iLon,iLat,-1:nAlts + 2,iIon,iBlock) =  1.0e6

              IDensityS(iLon,iLat,-1:nAlts + 2,ie_,iBlock) = &
                   IDensityS(iLon,iLat,-1:nAlts + 2,ie_,iBlock) + &
                   IDensityS(iLon,iLat,-1:nAlts + 2,iIon,iBlock) 

           enddo ! end inner ispecies loop

           !
           ! End ion loop
           !/


        enddo! end iLon loop
     enddo ! end iLat loop


     write(*,*) '============> init_msis.Mars.f90 Major Diagnostics:  Begin'

     where(NDensityS < 1.0e+03)
        NDensityS = 1.0e+03
     end where

     where(IDensityS < 1.0e+03)
        IDensityS = 1.0e+03
     end where

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


! do iSpecies = 1, nSpeciesTotal
!   write(*,*) 'iSpecies =', iSpecies
!
!   do iAlt = -1, nAlts + 2
!     write(*,*) 'iAlt =', iAlt
!
!     do iLat = -1, nLats + 2
!       write(*,*) 'iLat =', iLat
!
!       do iLon = -1, nLons + 2
!         write(*,*) 'iLon =', iLon
!
!        write(*,*) 'NDensityS(', iLon, iLat, iAlt, iSpecies,') =', &
!                    NDensityS(iLon,iLat,iAlt,iSpecies,iBlock)
!       write(*,*) 'Temperature(', iLon, iLat, iAlt,') =', &
!                   Temperature(iLon,iLat,iAlt,iBlock)
!
!       enddo
!     enddo
!   enddo
! enddo


!\
! Diagnostic Outputs
!/
!112 FORMAT(F6.1,1X,F8.4,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3, &
!              1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3, &
!              1X,ES10.3,1X,ES10.3)
!
!  open(UNIT = 26, FILE = 'test_neutrals.txt', STATUS='NEW',ACTION = 'WRITE')
!
!  do iiAlt = -1,nAlts + 2 
!        write(26,112) &
!      	newalt(iiAlt), &
!       	Temperature(iiLon,iiLat,iiAlt,iBlock), &
!       	NDensityS(iiLon,iiLat,iiAlt,iH_,iBlock), &
!       	NDensityS(iiLon,iiLat,iiAlt,iH2_,iBlock), &
!       	NDensityS(iiLon,iiLat,iiAlt,iCH_,iBlock), &
!       	NDensityS(iiLon,iiLat,iiAlt,i1CH2_,iBlock), &
!       	NDensityS(iiLon,iiLat,iiAlt,i3CH2_,iBlock), &
!       	NDensityS(iiLon,iiLat,iiAlt,iCH3_,iBlock), &
!       	NDensityS(iiLon,iiLat,iiAlt,iCH4_,iBlock), &
!       	NDensityS(iiLon,iiLat,iiAlt,iC2H4_,iBlock), &
!       	NDensityS(iiLon,iiLat,iiAlt,iN4S_,iBlock), &
!       	NDensityS(iiLon,iiLat,iiAlt,iN2_,iBlock), &
!       	NDensityS(iiLon,iiLat,iiAlt,iHCN_,iBlock), &
!       	NDensityS(iiLon,iiLat,iiAlt,iH2CN_,iBlock), &
!  end do
!
!  close(Unit = 26)


!\
! Calculating MeanMajorMass -----------------------------
!/

!\
! Initialize MeanMajorMass to 0.0
!/

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

     write(*,*) '==> Now Completing Mars Background Composition: END', iBlock   

  enddo

end subroutine init_msis

subroutine msis_bcs(iJulianDay,UTime,Alt,Lat,Lon,Lst, &
             F107A,F107,AP,LogNS, Temp, LogRho)

  write(*,*) "You can not use MSIS with any planet except Earth!!!"
  write(*,*) "If you ARE running Earth, then make the code again, using"
  write(*,*) "configure Earth ; make"
  call stop_gitm("I can not continue...")

end subroutine msis_bcs
