!Just compute EUVAC fluxes (12/16/2020)
subroutine bponder_euvac
  use ModConstants
  use ModInputs
  use ModIoUnit, only : UnitTmp_
  use ModEUV
  implicit none
  
  real :: P 
  real :: f107_bp !W/Hz/m2
  real :: f107_ave_bp

  logical :: exist
  integer :: iWave

  f107_bp = f107
  f107_ave_bp = f107
 
  P = (f107_bp + f107_ave_bp)/2.0 
  
  photonFlux_bp = f74113_bp*(1.0 + afac_bp*(P - 80))
  !f74113 = photons/m2/sec
  !P = F10.7 = W * sec/m^2

  !final units = photons * W / m^4
  
  do iWave = 1,nWavelengths
    if (photonFlux_bp(iWave) < 0.0) photonFlux_bp(iWave) = 0.0
  enddo

  photonFlux_bp = photonFlux_bp/(SunPlanetDistance**2) 
  
  
end subroutine bponder_euvac


!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine euv_ionization_heat(iBlock)

  use ModGITM
  use ModEUV
  use ModPlanet
  use ModConstants
  use ModInputs
  use ModSources
  use ModTime, only : tSimulation, CurrentTime, iTimeArray, iStep
  use ModIndicesInterfaces
  
  implicit none

  integer, intent(in) :: iBlock

  integer :: iAlt, iWave, iSpecies, iNeutral, iIon, iError, iLon,iLat
  integer :: iRxn
  real, dimension(nLons,nLats) :: Tau, Intensity, Intensity_bp

  logical :: IsFirstTime(nBlocksMax) = .true.

  real :: NeutralDensity(nLons, nLats, nSpeciesTotal)
  real :: ChapmanLittle(nLons, nLats, nSpecies)
  real :: EHeat(nLons, nLats), EHeat_bp(nLons, nLats), IHeat(nLons,nLats)
  real :: nEuvHeating(nLons,nLats,nAlts), neEuvHeating(nLons,nLats,nAlts)
  
  real :: wavelength_ave
  logical :: exist
  real :: TauS(nWavelengths,nAlts)
  !real :: TauTotal(nAlts)

  if (IsFirstTime(iBlock)) then

     IsFirstTime(iBlock) = .false.

     ! This transfers the specific photo absorption and ionization cross
     ! sections into general variables, so we can use loops...

     call fill_photo
     
     call init_energy_required_for_dissociation
     !if (iProc .eq. 1) then
     !  do iSpecies = 1,nSpeciesTotal
     !     write(*,*) cSpecies(iSpecies)
     !  enddo 
     !  write(*,*) energyRequired/1.602e-19
     !
     !endif
     !Manually hack in the CO2 dissociation threshold, because it's a bit smaller
     !Shaw et al., 1995
     energyRequired(iCO2_) = 5.453*1.602e-19 !eV -> correct units
 
  else
     if (floor((tSimulation - dT)/dTAurora) == &
          floor(tSimulation/dTAurora)) return
  endif

  call report("euv_ionization_heat",2)
  
  call start_timing("euv_ionization_heat")

  ChapmanLittle = 0.0
  EuvIonRate = 0.0
  EuvHeating(:,:,:,iBlock)= 0.0
  EuvHeating_bp(:,:,:,iBlock) = 0.0
  PhotoElectronHeating(:,:,:,iBlock)= 0.0
  eEuvHeating(:,:,:,iBlock) = 0.0
  EuvIonRateS(:,:,:,:,iBlock) = 0.0
  EuvDissRateS(:,:,:,:,iBlock) = 0.0
  EuvDissRateCO2_at_wave(:,:,:,:,iBlock) = 0.0
  nEuvHeating(:,:,:) = 0.0
  neEuvHeating(:,:,:) = 0.0
  DissociationHeatingRate(:,:,:,iBlock) = 0.0

  call bponder_euvac

  !flux_of_euv = photonFlux_bp
  call chapman_integrals(iBlock)
  
  do iAlt = 1, nAlts

     NeutralDensity = NDensityS(1:nLons,1:nLats,iAlt,1:nSpeciesTotal,iBlock)
     ChapmanLittle  = Chapman(:,:,iAlt,1:nSpecies,iBlock)
     EHeat = 0.0     
     EHeat_bp = 0.0
     IHeat = 0.0
     do iWave = 1, nWavelengths
        Tau = 0.0
        TauS(iWave,iAlt) = 0.0
       
        do iSpecies = 1, nSpecies
           Tau = Tau + &
                PhotoabsorptionCrossSection(iWave, iSpecies) * &
                ChapmanLittle(:,:,iSpecies)

           !Save the noon taus
           if (iSpecies .eq. nSpecies .and. iProc .eq. 2) then 
              TauS(iWave,iAlt) = Tau(3,1)

              !if (iWave .eq. nWavelengths) then
              !   TauTotal(iAlt) = SUM(TauS(:,iAlt))
              !endif
           endif
        enddo

        
        
        Intensity = Flux_of_EUV(iWave) * exp(-1.0*Tau)
        Intensity_bp = photonFlux_bp(iWave)*exp(-1.0*Tau)


        do iIon = 1, nIons-1
           iNeutral = PhotoIonFrom(iIon)
           EuvIonRateS(:,:,iAlt,iIon,iBlock) = &
                EuvIonRateS(:,:,iAlt,iIon,iBlock) + &
                Intensity * photoionizationCrossSection(iWave,iIon) * &
                NeutralDensity(:,:,iNeutral) * &
                (1.0 + PhotoElecIon(iWave,iIon))

        enddo

        wavelength_ave = (shortWavelengths(iWave) + longWavelengths(iWave))/2.0
        PhotonEnergy(iWave)= Planck_Constant*Speed_Light/wavelength_ave
        
        do iSpecies = 1, nSpecies !was nSpeciesTotal, but removing for photoDis
           EuvDissRateS(:,:,iAlt,iSpecies,iBlock) = &
                EuvDissRateS(:,:,iAlt,iSpecies,iBlock) + &
                Intensity*photoDissociationCrossSection(iWave,iSpecies) * &
                NeutralDensity(:,:,iSpecies) * &
                (1.0 + PhotoElecDiss(iWave,iSpecies))
           
           !Excess dissociation heating
           !Only applicable for CO2 for now
           !Starts when the energy > dissociation threshold, stops when
           !energies start entering into ionization-size levels
           !I believe EuvDissRateS does not work for this and 
           !a specific dissociation production rate AT A SPECIFIC wavelength
           !is needed instead. Hence the different variable.
           !Rho*cp adjustment below...           
           if (PhotonEnergy(iWave) - energyRequired(iSpecies) > 0.0 &
               .and. iSpecies == iCO2_ .and. useDissociationHeating) then
              if (photoionizationCrossSection(iWave,iCO2P_) .eq. 0.0) then
                 EuvDissRateCO2_at_wave(1:nLons,1:nLats,iAlt,iSpecies,iBlock) = &
                      Intensity*photoDissociationCrossSection(iWave,iSpecies) * &
                      NeutralDensity(:,:,iSpecies)
                 
                 DissociationHeatingRate(1:nLons,1:nLats,iAlt,iBlock) = &
                   DissociationHeatingRate(1:nLons,1:nLats,iAlt,iBlock) + &
                   2.0/3.0 * &
                   (PhotonEnergy(iWave) - energyRequired(iSpecies)) * &
                   EuvDissRateCO2_at_wave(1:nLons,1:nLats,iAlt,iSpecies,iBlock)
              endif
           endif
        enddo

        do iSpecies = 1, nSpecies
          if (shortWavelengths(iWave) < 1.0e-6) then
           EHeat = EHeat + &
                Intensity*PhotonEnergy(iWave)* &
                PhotoabsorptionCrossSection(iWave, iSpecies) * NeutralDensity(:,:,iSpecies)

           EHeat_bp = EHeat_bp + &
                Intensity_bp*PhotonEnergy(iWave)* &
                PhotoabsorptionCrossSection(iWave, iSpecies) * NeutralDensity(:,:,iSpecies)
           else if (shortWavelengths(iWave) .ge. 1.0e-6) then !this is 2.7 and 4.3
             IHeat = IHeat + &
                Intensity*PhotonEnergy(iWave)* &
                PhotoabsorptionCrossSection(iWave, iSpecies) * NeutralDensity(:,:,iSpecies)
           endif
        enddo
     enddo
     
     EuvHeating(:,:,iAlt,iBlock)  = EHeat*HeatingEfficiency_CB(:,:,iAlt,iBlock)
     EuvHeating_bp(:,:,iAlt,iBlock) = EHeat_bp*HeatingEfficiency_CB(:,:,iAlt,iBlock)
     eEuvHeating(:,:,iAlt,iBlock) = EHeat*eHeatingEfficiency_CB(:,:,iAlt,iBlock)
     
     if (useIRHeating) then
       QnirTOT(:,:,iAlt,iBlock) = IHeat*HeatingEfficiency_IR(:,:,iAlt,iBlock)
     else
       QnirTOT(:,:,iAlt,iBlock) = 0.0     
     endif

     do ilon = 1, nlons 
        do ilat =1 ,nlats
           if (Altitude_GB(iLon,iLat,iAlt,iBlock) .lt. 80000.0) then
              EUVHeating(iLon,iLat,iAlt,iBlock) =0.0
              eEUVHeating(iLon,iLat,iAlt,iBlock) =0.0
              EUVHeating_BP(iLon,iLat,iAlt,iBlock) = 0.0
           endif
        enddo
     enddo
  enddo
  

  !debugging
  !if (iStep > 500) then
  !  do iLon = 1,nLons
  !    do iLat = 1,nLats
  !      !write(*,*)                                                                                   
  !      if (sza(iLon,iLat,iBlock)*180.0/pi > 89.0 .and. &
  !          sza(iLon,iLat,iBlock)*180.0/pi < 95.0 .and. &
  !          latitude(iLat,iBlock) * 180.0/pi > 0.0 .and. &
  !          latitude(iLat,iBlock) * 180.0/pi < 1.5) then
  !         write(*,*) "SZA:", SZA(iLon,iLat,iBlock)*180.0/pi
  !         write(*,*) "Longitude:", longitude(iLon,iBlock)*180.0/pi
  !         write(*,*) "QnirTOT", QnirTOT(iLon,iLat,30:60,iBlock)*86400.0
  !         write(*,*) "EUV", EuvHeating(iLon,iLat,90,iBlock)*86400.0 / &
  !                           Rho(iLon,iLat,90,iBlock) / &
  !                           cp(iLon,iLat,90,iBlock)
  !       endif
  !    enddo
  !  enddo
  !  call stop_gitm("Debugging in calc_euv.f90")
  !endif


  if (IncludeEclipse) call calc_eclipse_effects
  if (IsEarth) call night_euv_ionization

  EuvIonRateS = EuvIonRateS + nEuvIonRateS
  EuvHeating(:,:,:,iBlock) = EuvHeating(:,:,1:nAlts,iBlock) + nEuvHeating
  EuvHeating_bp(:,:,:,iBlock) = EuvHeating_bp(:,:,1:nAlts,iBlock) + nEuvHeating
  eEuvHeating(:,:,:,iBlock) = eEuvHeating(:,:,1:nAlts,iBlock) + neEuvHeating

  !\
  ! Zero out EuvHeating if specified not to use it.
  !/

  if (UseSolarHeating) then

     EuvHeating2d = 0.0

     !EuvHeating = EuvHeating_bp
     do iAlt = 1, nAlts

        EuvHeating2d(1:nLons,1:nLats) = &
             EuvHeating2d(1:nLons,1:nLats) + &
             EuvHeating(:,:,iAlt,iBlock)  * &
             dAlt_GB(1:nLons,1:nLats,iAlt,iBlock)

        EuvHeating(1:nLons,1:nLats,iAlt,iBlock) = EuvHeating(1:nLons,1:nLats,iAlt,iBlock) / &
             Rho(1:nLons,1:nLats,iAlt,iBlock) / &
             cp(1:nLons,1:nLats,iAlt,iBlock) / &
             TempUnit(1:nLons,1:nLats,iAlt)


        EuvHeating_bp(1:nLons,1:nLats,iAlt,iBlock) = &
             EuvHeating_bp(1:nLons,1:nLats,iAlt,iBlock) / &
             Rho(1:nLons,1:nLats,iAlt,iBlock) / &
             cp(1:nLons,1:nLats,iAlt,iBlock) / &
             TempUnit(1:nLons,1:nLats,iAlt)

        DissociationHeatingRate(1:nLons,1:nLats,iAlt,iBlock) = &
             DissociationHeatingRate(1:nLons,1:nLats,iAlt,iBlock) / &
             Rho(1:nLons,1:nLats,iAlt,iBlock) / &
             cp(1:nLons,1:nLats,iAlt,iBlock) / &
             TempUnit(1:nLons,1:nLats,iAlt)
        
        EuvTotal(:,:,iAlt,iBlock) = EuvHeating(:,:,iAlt,iBlock) * &
             TempUnit(1:nLons,1:nLats,iAlt) / &
             HeatingEfficiency_CB(:,:,iAlt,iBlock)

        QnirTOT(1:nLons,1:nLats,iAlt,iBlock) = QnirTOT(1:nLons,1:nLats,iAlt,iBlock) / &
             Rho(1:nLons,1:nLats,iAlt,iBlock) / &
             cp(1:nLons,1:nLats,iAlt,iBlock) / &
             TempUnit(1:nLons,1:nLats,iAlt)

     enddo
  else
     EuvHeating = 0.0
  endif

  do iLon = 1,nLons
    do iLat = 1, nLats
      do iAlt = 1,nAlts
        if (DissociationHeatingRate(iLon,iLat,iAlt,iBlock) < -1.0) then
          !write(*,*) "euvDiss", EuvDissRateS(iLon,iLat,iAlt,:,iBlock)
          !write(*,*) "Dissociation value:", DissociationHeatingRate(iLon,iLat,iAlt,iBlock)
          !write(*,*) "iProc, iLon, iLat, iAlt", iProc, iLon,iLat,iAlt
          DissociationHeatingRate(iLon,iLat,iAlt,iBlock) = 0.0
        endif
       
      end do
    enddo
  enddo

  !Find 12 LST lon and lat
  !do iLon = 1,nLons
  !   do iLat = 1,nLats
  !     !write(*,*) longitude(iLon,iBlock)*180.0/pi
  !     if (longitude(iLon,iBlock)*180.0/pi > 70.0 .and. &
  !       longitude(iLon,iBlock)*180.0/pi < 83.0 .and. &
  !       latitude(iLat,iBlock)*180.0/pi < 11.0 .and. &
  !       latitude(iLat,iBlock)*180.0/pi > -10.0) then
  !        write(*,*) "iLon, iLat, iProc", iLon, iLat, iProc
  !        write(*,*) longitude(iLon,iBlock)*180.0/pi, latitude(iLat,iBlock)*180.0/pi
  !      endif
  !    enddo
  ! enddo
  
  
  !call stop_gitm("Found noon in calc_euv.f90")
  
  !Write Tau to a file
  exist = .false.
  if (iProc == 2 .and. .True.) then  
       inquire(file="taus.txt", exist=exist)
       if (exist) then
         open(84, file="taus.txt", status="old", position="append", action="write")
       else
         open(84, file="taus.txt", status="new", action="write")
      end if

      !Write a header with the columns
      ! Alt (km) Wvlgth 1 Wvlgth 2 ... sum(wavelengths)
      !wvavg(iWave)*1.0e10
      
      write(*,*) (shortWavelengths(1) + longWavelengths(1))*1.0e10/2
      write(*,*) (shortWavelengths(10) + longWavelengths(10))*1.0e10/2
      write(*,*) (shortWavelengths(20) + longWavelengths(20))*1.0e10/2
      write(*,*) (shortWavelengths(30) + longWavelengths(30))*1.0e10/2
      write(*,*) (shortWavelengths(41) + longWavelengths(41))*1.0e10/2
      
      
      write(84, *) &
           "Alt. (km), 1.5, 175.0, 425.0, 675.0, 975"
           
      
       do iAlt = 1,nAlts
         write(84, "(F5.1, 1X," // &
              "ES9.3, 1X," // &
              "ES9.3, 1X," // &
              "ES9.3, 1X," // &
              "ES9.3, 1X," // &
              "ES9.3)") &
         Altitude_GB(3,1,iAlt,iBlock)/1000.0, &
         TauS(1,iAlt), &
         TauS(10,iAlt), &
         TauS(20,iAlt), &
         TauS(30,iAlt), &
         TauS(41,iAlt)
         !TauTotal(iAlt)
       enddo
       close(84)
       call stop_gitm("In Calc_euv.f90")
   endif
  
  
  call end_timing("euv_ionization_heat")

contains 

  !\
  ! -------------------------------------------------------------
  !/

  subroutine calc_eclipse_effects

    real :: x,y,z,xp,yp,zp, distance, PercentDone
    real :: yPos, zPos, e
    real :: factor(nLons,nLats,nAlts)

    Factor = 1.0

    if ( CurrentTime > EclipseStartTime .and. &
         CurrentTime < EclipseEndTime) then

       PercentDone = (CurrentTime - EclipseStartTime) / &
            (EclipseEndTime - EclipseStartTime)

       yPos = EclipseStartY + (EclipseEndY - EclipseStartY) * PercentDone
       zPos = EclipseStartZ + (EclipseEndZ - EclipseStartZ) * PercentDone

       ! Need to rotate lat/lon/alt into solar coordinate X,Y,Z system

       !   x = r * cos(localtime-!pi) * cos(lats*!dtor)
       !   y = r * sin(localtime-!pi) * cos(lats*!dtor)
       !   z = r * sin(lats*!dtor)
       !
       !   xp = x * cos(tilt) - z * sin(tilt)
       !   zp = x * sin(tilt) + z * cos(tilt)

       do iLon = 1,nLons 
          do iLat = 1,nLats
             do iAlt = 1,nAlts
                ! Calculate X,Y,Z for points on Earth (in lat/localtime/alt)
                x = &
                     RadialDistance_GB(iLon,iLat,iAlt,iBlock) * &
                     cos(localtime(iLon)*pi/12.0 - pi) * &
                     cos(Latitude(iLat,iBlock))
                y = &
                     RadialDistance_GB(iLon,iLat,iAlt,iBlock) * &
                     sin(localtime(iLon)*pi/12.0 - pi) * &
                     cos(Latitude(iLat,iBlock))
                z = &
                     RadialDistance_GB(iLon,iLat,iAlt,iBlock) * &
                     sin(Latitude(iLat,iBlock))

                ! rotate around y-axis for the Earth Tilt:
                xp = x * cos(-SunDeclination) - z * sin(-SunDeclination)
                yp = y
                zp = x * sin(-SunDeclination) + z * cos(-SunDeclination)

                distance = sqrt( (yPos - yp)**2 + (zPos-zp)**2 )

                ! If it is close to the eclipse, do calculations

                if (Distance < EclipseMaxDistance + EclipseExpWidth*2) then

                   ! Assume a linearly decreasing effect:
                   e = EclipsePeak - EclipsePeak * distance/EclipseMaxDistance
                   if (e < 0) e = 0

                   ! Towards the edge, assume an exponential:
                   e = e + EclipseExpAmp * exp(-((Distance-EclipseMaxDistance)/EclipseExpWidth)**2)
                   Factor(iLon,iLat,iAlt) = 1.0 - e

                   ! Change the ionization rate and heating rate by the factor:
                   do iIon = 1, nIons-1
                      EuvIonRateS(iLon,iLat,iAlt,iIon,iBlock) = &
                           EuvIonRateS(iLon,iLat,iAlt,iIon,iBlock) * &
                           Factor(iLon,iLat,iAlt)
                   enddo
                   EuvHeating(iLon,iLat,iAlt,iBlock) = &
                        EuvHeating(iLon,iLat,iAlt,iBlock) * &
                        Factor(iLon,iLat,iAlt)

                   N2PERateS(iLon,iLat,iAlt,:,iBlock) = &
                        N2PERateS(iLon,iLat,iAlt,:,iBlock) * Factor(iLon,iLat,iAlt)
                   O2PERateS(iLon,iLat,iAlt,:,iBlock) = &
                        O2PERateS(iLon,iLat,iAlt,:,iBlock) * Factor(iLon,iLat,iAlt)
                   OPERateS(iLon,iLat,iAlt,:,iBlock) = &
                        OPERateS(iLon,iLat,iAlt,:,iBlock) * Factor(iLon,iLat,iAlt)
                   CH4PERateS(iLon,iLat,iAlt,:,iBlock) = &
                        CH4PERateS(iLon,iLat,iAlt,:,iBlock) * Factor(iLon,iLat,iAlt)

                endif

             enddo
          enddo
       enddo

    endif

  end subroutine calc_eclipse_effects

  !\
  ! -------------------------------------------------------------
  !/

  subroutine night_euv_ionization

!!! Nighttime EUV Ionization

    real :: night_col(nLons,nLats,nAlts,nSpecies)
    real :: dAlts(nLons,nLats)
    real, dimension(nLons,nLats) :: nTau, nIntensity
    real :: nNeutralDensity(nLons, nLats, nSpecies)
    real :: tempflux(nLons,nLats), SZALocal(nLons,nLats)
    real :: ScaleHeightS(nLons,nLats)
    real :: scatterfac = 0.01
    integer :: sinexpo = 15

    real :: nEHeat(nLons, nLats)
    real :: nPhotonEnergy(Num_NightWaveLens)
    real :: tempPhotonEnergy = 0.0

    nEuvIonRates   = 0.0
    night_col   = 0

    SZALocal = SZA(1:nLons,1:nLats,iBlock)

    do iWave=1, Num_NightWaveLens
       tempflux = 0.0
       select case (iWave)
       case (1)
          where(SZALocal .GE. Pi/2.) tempflux = -0.573563 * SZALocal**3 &
               + 4.66159 * SZALocal**2 - 12.6936 * SZALocal +11.66

          where(SZALocal .LT. Pi/2.) tempflux = (-0.573563 * (Pi/2.)**3 &
               + 4.66159 * (Pi/2.)**2 - 12.6936 * (Pi/2.) +11.66) &
               * sin(SZALocal)**sinexpo

          tempflux = tempflux * Flux_of_EUV(18) * scatterfac

          where(tempflux .LE. 0.) tempflux = 0.0

          tempPhotonEnergy = PhotonEnergy(18)

       case (2)
          where(SZALocal .GE. Pi/2.) tempflux = - 4.39459 * SZALocal**3 &
               + 28.8112*SZALocal**2  - 62.9068*SZALocal + 45.7575

          where(SZALocal .LT. Pi/2.) tempflux = (- 4.39459 * (Pi/2.)**3 &
               + 28.8112*(Pi/2.)**2  - 62.9068*(Pi/2.)  + 45.7575) &
               * sin(SZALocal)**sinexpo         

          tempflux = tempflux * Flux_of_EUV(35) * scatterfac

          where(tempflux .LE. 0.0) tempflux = 0.0

          tempPhotonEnergy = PhotonEnergy(35)

       case (4)    
          where(SZALocal .GE. Pi/2.) tempflux = - 0.356061 * SZALocal**3 &
               + 3.0883 * SZALocal**2  - 8.98066 *SZALocal+ 8.86673

          where(SZALocal .LT. Pi/2.) tempflux = (- 0.356061 * (Pi/2.)**3 &
               + 3.0883 * (Pi/2.)**2  - 8.98066 *(Pi/2.) + 8.86673) &
               * sin(SZALocal)**sinexpo

          tempflux = tempflux * Flux_of_EUV(12) * scatterfac

          where(tempflux .LE. 0.0) tempflux = 0.0

          tempPhotonEnergy= PhotonEnergy(12) 

       case default
          tempflux = 0.0
          tempPhotonEnergy = 0.0
       end select

       nighteuvflux(iWave,:,:,iBlock) = tempflux
       nPhotonEnergy(iWave) = tempPhotonEnergy

    end do


    !! Calculate Tau
    do iAlt= nAlts, 1, -1
       nNeutralDensity = NDensityS(1:nLons,1:nLats,iAlt,1:nSpecies,iBlock)
       dAlts = dAlt_GB(1:nLons,1:nLats,iAlt,iBlock)

       do iSpecies = 1, nSpecies

          if (iAlt .LT. nAlts) then

             night_col(:,:,iAlt,iSpecies) = night_col(:,:,iAlt+1,iSpecies) &
                  + nNeutralDensity(:,:,iSpecies)*dAlts

          else

             ScaleHeightS =  &
                  Temperature(1:nLons,1:nLats,iAlt,iBlock) &
                  * TempUnit(1:nLons,1:nLats,iAlt) * Boltzmanns_Constant &
                  / (-Gravity_GB(1:nLons,1:nLats,iAlt,iBlock) * Mass(iSpecies))
             night_col(:,:,iAlt,iSpecies) = nNeutralDensity(:,:,iSpecies) &
                  * ScaleHeightS
          endif

       enddo
    enddo


    do iAlt = 1, nAlts
       nEHeat = 0.0
       nNeutralDensity = NDensityS(1:nLons,1:nLats,iAlt,1:nSpecies,iBlock)
       do iWave = 1, Num_NightWaveLens
          nTau = 0.0

          do iSpecies = 1, nSpecies
             nTau = nTau + &
                  night_photoabs(iWave,iSpecies) * night_col(:,:,iAlt,iSpecies)
          enddo

          nIntensity = nighteuvflux(iWave,:,:,iBlock) * exp(-1.0*nTau)

          do iIon = 1, nIons-1
             iNeutral = PhotoIonFrom(iIon)
             nEuvIonRateS(:,:,iAlt,iIon,iBlock) = &
                  nEuvIonRateS(:,:,iAlt,iIon,iBlock) + &
                  nIntensity*night_photoion(iWave,iIon) * &
                  nNeutralDensity(:,:,iNeutral)
          enddo

          do iSpecies = 1, nSpecies
             nEHeat = nEHeat + &
                  nIntensity*nPhotonEnergy(iWave)* &
                  night_photoabs(iWave, iSpecies) * nNeutralDensity(:,:,iSpecies)            

          enddo

       enddo

       nEuvHeating(:,:,iAlt)  = nEHeat*HeatingEfficiency_CB(:,:,iAlt,iBlock)
       neEuvHeating(:,:,iAlt) = nEHeat*eHeatingEfficiency_CB(:,:,iAlt,iBlock)
       do ilon = 1, nlons 
          do ilat =1 ,nlats
             if (Altitude_GB(iLon,iLat,iAlt,iBlock) .lt. 80000.0) then
                nEUVHeating(iLon,iLat,iAlt) =0.0
                neEUVHeating(iLon,iLat,iAlt) =0.0
             endif
          enddo
       enddo

    enddo

  end subroutine night_euv_ionization



end subroutine euv_ionization_heat


!-------------------------------------------------------------------
! Subroutine for calculating the EUV flux in a vacuum.
!-------------------------------------------------------------------

subroutine calc_euv

  use ModEUV
  use ModInputs
 

  implicit none

  integer :: i, ii
  real    :: flxfac, wavelength_ave

  !:::::::::::::::::::::::::::::::: EUVAC :::::::::::::::::::::::
  !------ This EUV flux model uses the F74113 solar reference spectrum and
  !------ ratios determined from Hinteregger's SERF1 model. It uses the daily
  !------ F10.7 flux (F107) and the 81 day mean (F107A) as a proxy for 
  !------ scaling
  !------ The fluxes are returned in EUVFLX and correspond to the 37
  !------ wavelength
  !------ bins of Torr et al. [1979] Geophys. Res. Lett. p771.
  !------ See Richards et al. [1994] J. Geophys. Res. p8981 for details.
  !
  !...... F107   = input daily 10.7 cm flux index. (e.g. 74)
  !...... F107A  = input 81 day ave. of daily F10.7 centered on current day
  !...... EUVFLX = output array for EUV flux in units of photons/cm2/sec.
  !
  !
  !----- loop through the wavelengths calculating the scaling factors and
  !----- the resulting solar flux.
  !----- The scaling factors are restricted to be greater than 0.8
  !

  do i = 1, Num_waveLengths_Low

     ! This has to be backwards, due to the way that the wavelengths are defined:
     ii = Num_waveLengths_Low - i + 1
     FLXFAC=(1.0 + AFAC(II) * (0.5*(F107+F107A) - 80.0))
     IF(FLXFAC.LT.0.8) FLXFAC=0.8
     EUV_Flux(i) = F74113(II)* FLXFAC * 1.0E9 * 10000.
     
  enddo

end subroutine calc_euv

!-------------------------------------------------------------------
! Subroutine for calculating scaled solar flux.
!-------------------------------------------------------------------

subroutine calc_scaled_euv

  use ModEUV
  use ModInputs
  use ModTime
  use ModIndicesInterfaces
  use ModGITM, only : dt, iProc
  implicit none

  integer, parameter :: Hinteregger_Contrast_Ratio  = 0
  integer, parameter :: Hinteregger_Linear_Interp   = 1
  integer, parameter :: Tobiska_EUV91               = 2
  integer, parameter :: WoodsAndRottman_10Nov88     = 3
  integer, parameter :: WoodsAndRottman_20Jun89     = 4

  integer :: N, NN, iMin(1)=0,iError
  real    :: f107_Ratio, r1, r2, hlybr, fexvir, hlya, heiew
  real    :: xuvfac, hlymod, heimod, xuvf, wavelength_ave
  real (Real8_) :: rtime
  integer, dimension(7) :: Time_Array
  
 !DAVES:
  !real :: wvavg(Num_WaveLengths_High),SeeTime(nSeeTimes),tDiff(nSeeTimes)
  real :: wvavg(nWavelengths),SeeTime(nSeeTimes),tDiff(nSeeTimes)
  real :: y1(Num_WaveLengths_High), y2(Num_WaveLengths_High), x1, x2, x
  real :: m(Num_WaveLengths_High), k(Num_WaveLengths_High)
  character (len=2) :: dday, dhour, dminute 
  character (len=7) :: dtime

  real :: e
  integer :: iWave
  logical :: exist
  ! regression coefficients which reduce to solar min. spectrum:
  ! for Hinteregger_Contrast_Ratio model:

  integer :: newN
  real, dimension(1:3) :: B1, B2
  !real, dimension(Num_WaveLengths_High) :: Timed_Flux
  real, dimension(nWavelengths) :: Timed_Flux
  data B1/1.0, 0.0138, 0.005/
  data B2/1.0, 0.59425, 0.3811/

  ! 'best fit' regression coefficients, commented out, for reference:
  !     DATA B1/1.31, 0.01106, 0.00492/, B2/-6.618, 0.66159, 0.38319/

  iError = 0
  if (useRidleyEUV) then
    call get_f107(CurrentTime, f107, iError)
    if (iError /= 0) then
      write(*,*) "Error in getting F107 value.  Is this set?"
      write(*,*) "Code : ",iError
      call stop_gitm("Stopping in euv_ionization_heat")
   endif

    call get_f107a(CurrentTime, f107a, iError)
    if (iError /= 0) then
       write(*,*) "Error in getting F107a value.  Is this set?"
       write(*,*) "Code : ",iError
       call stop_gitm("Stopping in euv_ionization_heat")
    endif
  endif

  if (UseEUVData) then
     call Set_Euv(iError, CurrentTime, EndTime)

     do N=nWavelengths,1,-1
        !the bins are backwards in FISM
        wvavg(N)=(shortWavelengths(N) + longWavelengths(N))/2.
        !write(*,*) N, shortWavelengths(N), longWavelengths(N), wvavg(N)
     enddo

     call start_timing("new_euv")
     SeeTime(:) = 0

     do N = 1, nSeeTimes
        SeeTime(N) = TimeSee(N)
     enddo

     tDiff = CurrentTime - SeeTime
    
     where (tDiff .lt. 0) tDiff = 1.e20
     iMin = minloc(tDiff)
     if (iMin(1) .eq. 0) iMin(1) = 1
     
     Timed_Flux(nWavelengths-num_wavelengths_high+1:nWavelengths) = SeeFlux(1:num_wavelengths_high,iMin(1))

     if (CurrentTime .ge. FlareTimes(iFlare) .and. CurrentTime-dt .le. FlareTimes(iFlare)) then

        FlareStartIndex = iMin(1)+1
        FlareEndIndex = iMin(1) + FlareLength
        Timed_Flux = SeeFlux(:,FlareStartIndex)
        FlareEndTime = SeeTime(FlareEndIndex)
        DuringFlare = .true.
        iFlare = iFlare + 1
       
     else
        if (DuringFlare) then
           if (Seetime(iMin(1)+1) .lt. FlareTimes(iFlare) .or. FlareTimes(iFlare) .eq. 0) then
              if (CurrentTime .lt. SeeTime(FlareStartIndex)) then
                 !We may not be to the point where the flare has begun in the SEE data yet...
                 Timed_Flux(nWavelengths-num_wavelengths_high+1:nWavelengths) = SeeFlux(:,FlareStartIndex)
              else
                 if (CurrentTime .le. FlareEndTime) then
                    !Exponentially interpolate between last seetime and next seetim 
                    !using y = kexp(-mx)
                   
                    y1 = SeeFlux(:,iMin(1))
                    y2 = SeeFlux(:,iMin(1)+1)
                    x1 = 0
                    x2 = SeeTime(iMin(1)+1) - SeeTime(iMin(1))
                    x = CurrentTime - SeeTime(iMin(1))
                   
                    m = ALOG(y2/y1)/(x1-x2)
                    k = y1*exp(m*x1)
                    Timed_Flux = k*exp(-1*m*x)
                 else
                    DuringFlare = .False.
                 end if
              end if
           end if
        end if
     end if

     !These are for the last three wavelengths in euv.csv.
     !Hard coded fluxes from FISM2 on LASP site.
     !Timed_Flux(1800-1850ish) = (LASP site value)*5 nm
     !CSV file reads our wavelength bins...
     !multiply by dWavelength which is 5 nm.

     !2009 values
     !Timed_Flux(5) = 6.9e-3 !W/m^2 (because of dWavelength!)
     !Timed_Flux(4) = 1.03e-2
     !Timed_Flux(3) = 1.46e-2
     
     Timed_Flux(5) = fism175180 !W/m^2 (because of dWavelength!)                   
     Timed_Flux(4) = fism180185 !180-185nm                                              
     Timed_Flux(3) = fism185190 !185-190nm    

     !IR wavelength: previously 80% and we were too hot by ~50 K at 105 km. 
     !               trying 64% to see how this goes. 4.3 micron

     if (UseIRHeating .and. (.not. UseGilli) .and. (.not. UseRoldan)) then
       !V-GITM I values:
       !Timed_Flux(2) = 26530*0.64 !2.7 micron
       !Timed_Flux(1) = 26530*0.4 !temporary adjust to reduce peak heating, 4.3 micron     
 
       !Debugging values:
       Timed_Flux(2) = 0.18
       Timed_Flux(1) = 4.4 

     else
       Timed_Flux(1) = 0.0
       Timed_Flux(2) = 0.0
     endif

     !2002 values
     !Timed_Flux(60) = 7.9e-3 !W/m^2
     !Timed_Flux(61) = 1.14e-2
     !Timed_Flux(62) = 1.6e-2
                                             
     do N=1,nWavelengths
        Flux_of_EUV(N) = Timed_Flux(nWavelengths - (N-1))*wvavg(N)/(6.626e-34*2.998e8) &
             /(SunPlanetDistance**2)
     enddo

     call end_timing("new_euv")

  else
     if (UseRidleyEUV) then
        !what's the size of Solar_Flux
        !RidleyPowers???? data and size?
        !RidleySlopes same 
        !RidleyIntercepts same
        do N=nWavelengths,1,-1
          !the bins are backwards in FISM                                                        
          wvavg(N)=(shortWavelengths(N) + longWavelengths(N))/2.
          !write(*,*) N, shortWavelengths(N), longWavelengths(N), wvavg(N)                       
        enddo

        !N = 1 handles wvavg of 0.1 - 0.2 nm
        do N = 1, Num_Wavelengths_High
           Solar_Flux(N) = &
                RidleySlopes(1,N) * (f107**RidleyPowers(1,N)) + &
                RidleySlopes(2,N) * (f107a**RidleyPowers(2,N)) + &
                RidleySlopes(3,N) * (f107a-f107) + &
                RidleyIntercepts(N)

           !flux_of_euv is forward...
           Flux_of_EUV(N) = Solar_Flux(nWavelengths - (N-1))*wvavg(N)/(6.626e-34*2.998e8) &
                /(SunPlanetDistance**2)

        enddo
        
        !To occupy IR and 175-190 wavelength bins
        Flux_of_EUV(3:5) = 0.0

        if (UseIRHeating .and. (.not. UseGilli) .and. (.not. UseRoldan)) then
          Flux_of_EUV(nWavelengths-1) = 26530*0.64 * 2.7e-6/(6.626e-34*2.998e8)/(SunPlanetDistance**2)
          !2.7 micron                                                   
          Flux_of_EUV(nWavelengths) = 26530*0.4 * 4.3e-6/(6.626e-34*2.998e8)/(SunPlanetDistance**2)
          !temporary adjust to reduce peak heating, 4.3 micron        
        else
          Flux_of_EUV(nWavelengths-1:nWavelengths) = 0.0
        endif

        !what order is the waves.dat data? does it need to be reversed?

     else

        ! This runs EUVAC and puts the data into EUV_Flux
        if (UseEUVAC) call calc_euv

        ! This stuff is for Tobiska:
        hlybr = 0.
        fexvir = 0.
        hlya = 3.E+11 + 0.4E+10 * (f107-70.)
        heiew = 0.

        xuvfac = 4.0 - f107_ratio
        if (xuvfac < 1.0) xuvfac = 1.0

        ! I think that this is the Hinteregger stuff, which goes from
        ! 2A - 1750A, so it is used to fill in the Above and Below spectrum
        f107_ratio = (f107-68.0) / (243.0-68.0)
        do N = 1, Num_WaveLengths_High
           Solar_Flux(N) = RFLUX(N) + (XFLUX(N)-RFLUX(N)) * f107_Ratio
        enddo

        ! This is a total hack, just comparing FISM to these fluxes:
        Solar_flux(Num_WaveLengths_High-4) = Solar_flux(Num_WaveLengths_High-4) * 3.0
        Solar_flux(Num_WaveLengths_High-3) = Solar_flux(Num_WaveLengths_High-3) * 3.0
        Solar_flux(Num_WaveLengths_High-2) = Solar_flux(Num_WaveLengths_High-2) * 5.0
        Solar_flux(Num_WaveLengths_High-1) = Solar_flux(Num_WaveLengths_High-1) * 100.0
        Solar_flux(Num_WaveLengths_High) = Solar_flux(Num_WaveLengths_High) * 800.0

        if (.not.UseAboveHigh) Solar_flux(56:Num_WaveLengths_High) = 0
        if (.not.UseBelowLow)  Solar_flux(1:Num_WaveLengths_Low) = 0
  
        if (HLYA > 0.001) then

           hlymod = hlya

        else

           if (heiew > 0.001) then
              hlymod = heiew * 3.77847e9 + 8.40317e10
           else
              hlymod = 8.70e8 * F107 + 1.90e11
           endif

        endif

        if (heiew > 0.001) then
           heimod = heiew * 3.77847e9 + 8.40317e10
        else
           heimod = hlymod
        endif

        ! hlymod is SME Lyman-alpha
        ! heimod is He I 10,830A, scaled to Lyman-alpha
        do N=16,55
           Solar_Flux(N) = TCHR0(N)        + &
                TCHR1(N)*hlymod + &
                TCHR2(N)*heimod + &
                TCOR0(N)        + &
                TCOR1(N)*f107   + &
                TCOR2(N)*f107A
        enddo

        !
        ! Substitute in H Lyman-alpha and XUVFAC if provided:
        !

        if (hlya > 0.001) Solar_Flux(12) = hlya / 1.E9
        if (xuvfac > 0.001) THEN
           xuvf = xuvfac
        else
           xuvf = 1.0
        endif

        !
        ! Convert from gigaphotons to photons, cm^-2 to m^-2, etc.:
        !
        
        do N=1,nWavelengths

           IF (Solar_Flux(N) < 0.0) Solar_Flux(N) = 0.0

           ! I don't know why we scale this...
           IF ((longWavelengths(N)*1.0e10 < 251.0) .AND. &
               (shortWavelengths(N)*1.0e10 > 15.0)) then
              Solar_Flux(N) = Solar_Flux(N)*xuvf
           endif

           !
           ! Convert to photons/m^2/s
           !

           Solar_Flux(N) = Solar_Flux(N) * 1.E9 * 10000.0

           !
           ! Calculate the energy in the bin:
           !
     
           !the timed_flux bins are backwards from the wave the 
           !CSV file reads our wavelength bins...
           wavelength_ave = (shortWavelengths(59 - (N-1)) + longWavelengths(59 - (N-1)))/2.0
           PhotonEnergy(N)= 6.626e-34*2.998e8/(wavelength_ave)

        enddo
        ! Solar_Flux has the Tobiska and Hinteregger fluxes in there already:
  
        do N = 1,nWavelengths
           Flux_of_EUV(N) = Solar_Flux(N)
        enddo

        if (UseEUVAC) then
           do N=1,Num_WaveLengths_Low
              NN = N+15
              if (UseTobiska) then
                 Flux_of_EUV(NN) = 0.5*(EUV_Flux(N)+Solar_Flux(NN))
              else
                 Flux_of_EUV(NN) = EUV_Flux(N)
              endif
           enddo
        endif
  
        ! Take into account the sun distance to the planet:
        Flux_of_EUV = Flux_of_EUV/(SunPlanetDistance**2)
     endif

  endif
  
  !do iWave = 1,nWavelengths
    !if (abs(WAVES(iWave) - shortWavelengths(nWavelengths - (iWave-1))*1.0e10) > &
    !    1.0e-5) then
    !  write(*,*) "Diff in wavelength bins:", &
    !             WAVEL(iWave) - longWavelengths(nWavelengths - (iWave-1))*1.0e10
    !  write(*,*) "Diff in wavelengths bins (short):", &
    !             WAVES(iWave) - shortWavelengths(nWavelengths - (iWave-1))*1.0e10
    !endif
  !enddo
end subroutine calc_scaled_euv

 !-------------------------------------------------------------------
! Subroutine for initializing photoabsorption, photoionization,
! and branching ratio quantities.
!-------------------------------------------------------------------

subroutine init_euv

  use ModEUV
  use ModInputs, only: UseRidleyEUV

  implicit none

  integer :: N, NN, iError

  call report("init_euv",2)

  EUVEFF = 0.05

  do N = 1,15

     RLMSRC(N) = RLMEUV(N)

  enddo


  DO N = 1, Num_WaveLengths_Low

     NN = N+15    ! 16:52

     BranchingRatio_OPlus2P(N) = 0.
     IF (N.GT.14) then
        BranchingRatio_OPlus2P(N) =                                   &
             1. - BranchingRatio_OPlus2D(N) - BranchingRatio_OPlus4S(N)
     endif

     PhotoAbs_O2(NN)      = Photoabsorption_O2(N)
     PhotoAbs_O(NN)       = Photoabsorption_O(N) 
     PhotoAbs_N2(NN)      = Photoabsorption_N2(N)
     PhotoAbs_CO2(NN)   = Photoabsorption_CO2(N)
     PhotoAbs_CO(NN)   = Photoabsorption_CO(N)

     PhotoIon_O2(NN)      = Photoionization_O2(N)
     PhotoIon_CO2(NN)   = Photoionization_CO2(N)
     PhotoIon_CO(NN)   = Photoionization_CO(N)
     PhotoIon_OPlus4S(NN) = Photoionization_O(N)*BranchingRatio_OPlus4S(N)
     PhotoIon_N2(NN)      = Photoionization_N2(N)
     PhotoIon_N(NN)       = Photoionization_N(N)
     PhotoIon_OPlus2D(NN) = Photoionization_O(N)*BranchingRatio_OPlus2D(N)
     PhotoIon_OPlus2P(NN) = Photoionization_O(N)*BranchingRatio_OPlus2P(N)

     BranchingRatio_N2(NN) = BranchingRatio_N2_to_NPlus(N)
     BranchingRatio_O2(NN) = BranchingRatio_O2_to_OPlus(N)

  enddo

  ! Convert everything from /cm2 to /m2

  PhotoAbs_O2      = PhotoAbs_O2  / 10000.0
  PhotoAbs_O       = PhotoAbs_O   / 10000.0
  PhotoAbs_N2      = PhotoAbs_N2  / 10000.0
  PhotoAbs_CH4     = PhotoAbs_CH4 / 10000.0
  PhotoAbs_CO2      = PhotoAbs_CO2/ 10000.0
  PhotoAbs_CO      = PhotoAbs_CO  / 10000.0
  PhotoAbs_He       = PhotoAbs_He / 10000.0

  PhotoIon_O2      = PhotoIon_O2        / 10000.0
  PhotoIon_CO2      = PhotoIon_CO2        / 10000.0
  PhotoIon_CO      = PhotoIon_CO        / 10000.0
  PhotoIon_OPlus4S = PhotoIon_OPlus4S   / 10000.0
  PhotoIon_N2      = PhotoIon_N2        / 10000.0
  PhotoIon_N       = PhotoIon_N         / 10000.0
  PhotoIon_OPlus2D = PhotoIon_OPlus2D   / 10000.0
  PhotoIon_OPlus2P = PhotoIon_OPlus2P   / 10000.0
  
  if (UseRidleyEUV) call read_euv_waves(iError)
     
end subroutine init_euv

subroutine Set_Euv(iError, StartTime, EndTime)

  use ModKind
  use ModEUV
  use ModInputs
  use ModTime, only: iTimeArray

  implicit none

  integer, intent(out)     :: iError
  real(Real8_), intent(in) :: EndTime, StartTime

  character(len=20)                       :: line, cline
  integer, dimension(7)                   :: TimeOfFlare,TimeArray
  real, dimension(6+Num_Wavelengths_High) :: temp

  real (Real8_) :: TimeDelay, BufferTime = 180.0

  logical :: NotDone
  logical :: timesMatch
  integer :: i, iline, ioerror, nline, nILine = 1, iErrortemp = 0
  character(len = 10), dimension(30) :: date  
  integer, dimension(30) :: fismYear = 0.0
  integer, dimension(30) :: fismMonth = 0.0
  integer, dimension(30) :: fismDay = 0.0
  real, dimension(30) :: fism175180array = 0.0
  real, dimension(30) :: fism180185array = 0.0
  real, dimension(30) :: fism185190array = 0.0

  iError = 0

  ! If we have been here before and we read the entire file, leave
  if (.not.ReReadEUVFile .and. nSeeTimes > 0) return

  ! This if statement makes sure that we have been here before and 
  ! want to be here again
  if (ReReadEUVFile) then
     ! If we still have a lot of data in memory, then don't bother
     ! reading more.
     if (StartTime + BufferTime < TimeSee(nSeeTimes)) return
  endif

  ! Assume that we can read the entire file
  ReReadEUVFile = .false.

  cline = ' '

  open(unit = iInputUnit_, file=cEUVFile, IOSTAT = iError)
  

  if (iError /= 0) then
     write(*,*) "Error in opening EUV file  Is this set?"
     write(*,*) "Code : ",iError,cEUVFile
     call stop_gitm("Stopping in calc_euv")
  endif
     
  iline = 1

  NotDone = .true.
  do while (NotDone)
     read(iInputUnit_,'(a)',iostat=iError) cLine
     if (cline(1:1) .eq. '#') then 

        ! Remove anything after a space or TAB

        i=index(cLine,' '); if(i>0)cLine(i:len(cLine))=' '
        i=index(cLine,char(9)); if(i>0)cLine(i:len(cLine))=' '
        
        select case (cLine)
           
        case("#FLARES")
           read(iInputUnit_,*) nFlares
           do i = 1, nFlares 
              read(iInputUnit_,*) TimeOfFlare(1:6)
              TimeOfFlare(7) = 0
              call time_int_to_real(TimeOfFlare,FlareTimes(i))
           enddo
           
        case('#START')  
           NotDone = .false.
           
        end select
     end if
     
     iline = iline + 1
     
  enddo

  TimeSee(:) = 0
  iline = 1
  read(iInputUnit_,*,iostat=iError) temp

  iErrortemp = 0
  do while (iErrortemp .eq. 0)
     
     TimeArray(1:6) = temp(1:6)
     TimeArray(7) = 0
     call time_int_to_real(TimeArray,TimeSee(iLine))
     SeeFlux(:,iline) = temp(7:6+Num_WaveLengths_High)

     if ( TimeSee(iLine) > StartTime-BufferTime .and. &
          TimeSee(iLine) < EndTime+BufferTime) &
          iline = iline + 1

     read(iInputUnit_,*,iostat=iErrortemp) temp

     if (iline > nSeeLinesMax-1) then
        iErrortemp = 1
        ReReadEUVFile = .true.
     endif

  enddo

  close(iInputUnit_)
  nSeeTimes = iline - 1

  if (nSeeTimes .gt. 3) iError = 0

  !Read another file with FISM data for bins 175-190nm.
  !open(unit = iInputUnit_, file=cFISM2File, IOSTAT = iError)
  !if (iError /= 0) then
  !   write(*,*) "Error in opening fism_extraBins_YYYYmmdd.txt file  Is this set?"
  !   write(*,*) "Code : ",iError,cFISM2File
  !   call stop_gitm("Stopping in calc_euv -> set_euv")
  !endif

  !four lines of header info
  !read(iInputUnit_, *)
  !read(iInputUnit_, *)
  !read(iInputUnit_, *)
  !read(iInputUnit_, *)
  
  !This is only if the run is 10 days. Needs to be 
  !re-written in a correct fashion.
  !do iline = 1,10
  !  read(iInputUnit_,*) date(iline), fism175180array(iLine), &
  !                      fism180185array(iLine), &
  !                      fism185190array(iLine)
  !
  !
  !  read(date(iLine)(1:4), '(i4)') fismYear(iLine)
  !  read(date(iLine)(6:7), '(i2)') fismMonth(iLine)
  !  read(date(iLine)(9:10), '(i2)') fismDay(iLine)
  !enddo

  !close(iInputUnit_)

  !timesMatch = .False. 
  !do i = 1,size(fismYear)
  !  if (iTimeArray(1) .eq. fismYear(i) .and. &
  !      iTimeArray(2) .eq. fismMonth(i) .and. &
  !      iTimeArray(3) .eq. fismDay(i)) then
  !    timesMatch = .True.
  !    
  !    fism185190 = fism185190array(i)
  !    fism180185 = fism180185array(i)
  !    fism175180 = fism175180array(i)
  !
  !    exit
  !  endif
  !enddo

  !if (.not. timesMatch) call stop_gitm("Check UA/DataIn/fism_extraBins_????????.txt.")

end subroutine Set_Euv

  
subroutine read_euv_waves(iError)

  use ModKind
  use ModEUV
  use ModInputs

  implicit none

  integer, intent(out)     :: iError
  character(len=80)        :: cline

  integer :: nFiles, nWaves, iWave, i 
  real :: tmp(6)
  
  call report("read_euv_waves",2)

  open(unit = iInputUnit_, file=cRidleyEUVFile, IOSTAT = iError)

  if (iError /= 0) then
     write(*,*) "Error in opening RidleyEUVFile  Is this set?"
     write(*,*) "Code : ",iError,cRidleyEUVFile
     call stop_gitm("Stopping in read_euv_waves")
  endif
     
  read(iInputUnit_,*,iostat=iError) nFiles
  
  do i=1,nFiles
     read(iInputUnit_,*,iostat=iError) cline
  enddo

  read(iInputUnit_,*,iostat=iError) nWaves

  if (nWaves /= Num_WaveLengths_High) then
     write(*,*) "Error in RidleyEUVFile - nWaves /= Num_wavelengths_high"
     write(*,*) "Values : ", nWaves, Num_WaveLengths_High
     call stop_gitm("Stopping in read_euv_waves")
  endif

  read(iInputUnit_,*,iostat=iError) cline
  
  do i=1,nWaves

     read(iInputUnit_,*,iostat=iError) iWave, tmp
     RidleySlopes(1,i) = tmp(1)
     RidleySlopes(2,i) = tmp(2)
     RidleySlopes(3,i) = tmp(3)
     RidleyIntercepts(i) = tmp(4)
     RidleyPowers(1,i) = tmp(5)
     RidleyPowers(2,i) = tmp(6)
     
  enddo

  close(iInputUnit_)
  
end subroutine read_euv_waves
