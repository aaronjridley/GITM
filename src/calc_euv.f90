!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine euv_ionization_heat(iBlock)

  use ModGITM
  use ModEUV
  use ModPlanet
  use ModConstants
  use ModInputs
  use ModSources
  use ModTime, only : tSimulation, CurrentTime
  use ModIndicesInterfaces

  implicit none

  integer, intent(in) :: iBlock

  integer :: iAlt, iWave, iSpecies, iIon, iError, iLon,iLat
  real, dimension(nLons,nLats) :: Tau, Intensity

  logical :: IsFirstTime(nBlocksMax) = .true.

  real :: photoion(Num_WaveLengths_High, nIons-1)
  real :: photoabs(Num_WaveLengths_High, nSpeciesTotal)
  real :: photodis(Num_WaveLengths_High, nSpeciesTotal)

  real :: NeutralDensity(nLons, nLats, nSpecies)
  real :: ChapmanLittle(nLons, nLats, nSpecies)
  real :: EHeat(nLons, nLats)
  real :: nEuvHeating(nLons,nLats,nAlts), neEuvHeating(nLons,nLats,nAlts)

  if (IsFirstTime(iBlock)) then
     IsFirstTime(iBlock) = .false.
  else
     if (floor((tSimulation - dT)/dTAurora) == &
          floor(tSimulation/dTAurora)) return
  endif

  call report("euv_ionization_heat",2)
  call start_timing("euv_ionization_heat")
  call chapman_integrals(iBlock)

  EuvIonRate = 0.0
  EuvHeating(:,:,:,iBlock)= 0.0
  PhotoElectronHeating(:,:,:,iBlock)= 0.0
  eEuvHeating(:,:,:,iBlock) = 0.0
  EuvIonRateS(:,:,:,:,iBlock) = 0.0 
  EuvDissRateS(:,:,:,:,iBlock) = 0.0
  nEuvHeating(:,:,:) = 0.0
  neEuvHeating(:,:,:) = 0.0

  photoion(1:Num_Wavelengths_High,1:nIons-1) = 0.0
  photoabs(1:Num_Wavelengths_High,1:nSpeciesTotal)= 0.0
  photodis(1:Num_Wavelengths_High,1:nSpeciesTotal) = 0.0

  ! This transfers the specific photo absorption and ionization cross
  ! sections into general variables, so we can use loops...

  call fill_photo(photoion, photoabs, photodis)

  do iAlt = 1, nAlts

     NeutralDensity = NDensityS(1:nLons,1:nLats,iAlt,1:nSpecies,iBlock)
     ChapmanLittle  = Chapman(:,:,iAlt,1:nSpecies,iBlock)
     EHeat = 0.0

     do iWave = 1, Num_WaveLengths_High

        Tau = 0.0
        do iSpecies = 1, nSpecies
           Tau = Tau + &
                photoabs(iWave, iSpecies) * ChapmanLittle(:,:,iSpecies)
        enddo

        Intensity = Flux_of_EUV(iWave) * exp(-1.0*Tau)

        do iIon = 1, nIons-1
           EuvIonRateS(:,:,iAlt,iIon,iBlock) = &
                EuvIonRateS(:,:,iAlt,iIon,iBlock) + &
                Intensity*PhotoIon(iWave,iIon)
        enddo
 
       do iSpecies = 1, nSpeciesTotal
           EuvDissRateS(:,:,iAlt,iSpecies,iBlock) = &
                EuvDissRateS(:,:,iAlt,iSpecies,iBlock) + &
                Intensity*PhotoDis(iWave,iSpecies)
        enddo

        do iSpecies = 1, nSpecies
           EHeat = EHeat + &
                Intensity*PhotonEnergy(iWave)* &
                photoabs(iWave, iSpecies) * NeutralDensity(:,:,iSpecies)

        enddo

     enddo
     
     EuvHeating(:,:,iAlt,iBlock)  = EHeat*HeatingEfficiency_CB(:,:,iAlt,iBlock)
     eEuvHeating(:,:,iAlt,iBlock) = EHeat*eHeatingEfficiency_CB(:,:,iAlt,iBlock)
     do ilon = 1, nlons 
        do ilat =1 ,nlats
           if (Altitude_GB(iLon,iLat,iAlt,iBlock) .lt. 80000.0) then
              EUVHeating(iLon,iLat,iAlt,iBlock) =0.0
              eEUVHeating(iLon,iLat,iAlt,iBlock) =0.0
           endif
        enddo
     enddo
  enddo


  call night_euv_ionization

  EuvIonRateS = EuvIonRateS + nEuvIonRateS
  EuvHeating(:,:,:,iBlock) = EuvHeating(:,:,1:nAlts,iBlock) + nEuvHeating
  eEuvHeating(:,:,:,iBlock) = eEuvHeating(:,:,1:nAlts,iBlock) + neEuvHeating

  !\
  ! Zero out EuvHeating if specified not to use it.
  !/

  if (UseSolarHeating) then
     do iAlt = 1, nAlts

        EuvHeating(:,:,iAlt,iBlock) = EuvHeating(:,:,iAlt,iBlock) / &
           Rho(1:nLons,1:nLats,iAlt,iBlock) / &
           cp(1:nLons,1:nLats,iAlt,iBlock) / &
           TempUnit(1:nLons,1:nLats,iAlt)

        EuvTotal(:,:,iAlt,iBlock) = EuvHeating(:,:,iAlt,iBlock) * &
             TempUnit(1:nLons,1:nLats,iAlt) / &
             HeatingEfficiency_CB(:,:,iAlt,iBlock)

     enddo
  else
     EuvHeating = 0.0
  endif

  call end_timing("euv_ionization_heat")

contains 

  subroutine night_euv_ionization

!!! Nighttime EUV Ionization

    real :: night_photoion(Num_NightWaveLens, nIons-1)
    real :: night_photoabs(Num_NightWaveLens, nSpecies)
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

    night_photoion = 0.0
    night_photoabs = 0.0
    nEuvIonRates   = 0.0
    night_col   = 0
    night_photoion(:,iN2P_)    = Night_PhotoIon_N2 *1.e-4
    night_photoion(:,iO2P_)    = Night_PhotoIon_O2 *1.e-4
    night_photoion(:,iNOP_)    = Night_PhotoIon_NO *1.e-4
    night_photoion(:,iO_4SP_)  = Night_PhotoIon_OPlus4S *1.e-4
    night_photoion(:,iO_2DP_)  = Night_PhotoIon_OPlus2D *1.e-4
    night_photoion(:,iO_2PP_)  = Night_PhotoIon_OPlus2P *1.e-4
    night_photoion(:,iNP_)     = Night_PhotoIon_N *1.e-4
 
    night_photoabs(:,iN2_)    = Night_PhotoAbs_N2 *1.e-4
    night_photoabs(:,iO2_)    = Night_PhotoAbs_O2 *1.e-4
    night_photoabs(:,iNO_)    = Night_PhotoAbs_NO *1.e-4
    night_photoabs(:,iO_3P_)  = Night_PhotoAbs_O *1.e-4
    night_photoabs(:,iN_4S_)  = Night_PhotoAbs_N *1.e-4


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
             nEuvIonRateS(:,:,iAlt,iIon,iBlock) = &
                  nEuvIonRateS(:,:,iAlt,iIon,iBlock) + &
                  nIntensity*night_photoion(iWave,iIon)
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

  integer :: i
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

     FLXFAC=(1.0 + AFAC(I) * (0.5*(F107+F107A) - 80.0))
     IF(FLXFAC.LT.0.8) FLXFAC=0.8
     EUV_Flux(i) = F74113(I) * FLXFAC * 1.0E9 * 10000.

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
  use ModGITM, only : dt
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
  real :: wvavg(Num_WaveLengths_High),SeeTime(nSeeTimes),tDiff(nSeeTimes)
  real :: y1(Num_WaveLengths_High), y2(Num_WaveLengths_High), x1, x2, x
  real :: m(Num_WaveLengths_High), k(Num_WaveLengths_High)
  character (len=2) :: dday, dhour, dminute 
  character (len=7) :: dtime

  ! regression coefficients which reduce to solar min. spectrum:
  ! for Hinteregger_Contrast_Ratio model:

  real, dimension(1:3) :: B1, B2
  real, dimension(Num_WaveLengths_High) :: Timed_Flux
  data B1/1.0, 0.0138, 0.005/
  data B2/1.0, 0.59425, 0.3811/

  ! 'best fit' regression coefficients, commented out, for reference:
  !     DATA B1/1.31, 0.01106, 0.00492/, B2/-6.618, 0.66159, 0.38319/

  iError = 0
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

  if (UseEUVData) call Set_Euv(iError, CurrentTime, EndTime)

  call calc_euv

  hlybr = 0.
  fexvir = 0.
  hlya = 3.E+11 + 0.4E+10 * (f107-70.)
  heiew = 0.

  f107_ratio = (f107-68.0) / (243.0-68.0)

  xuvfac = 4.0 - f107_ratio
  if (xuvfac < 1.0) xuvfac = 1.0

  do N = 1, Num_WaveLengths_High
     Solar_Flux(N) = RFLUX(N) + (XFLUX(N)-RFLUX(N)) * f107_Ratio
  enddo

  iModelSolar = Tobiska_EUV91
!  iModelSolar = Hinteregger_Contrast_Ratio

  select case(iModelSolar)

  case (Hinteregger_Contrast_Ratio)

     if (hlybr > 0.001) then
        r1 = hlybr
     else
        r1 =  B1(1) + B1(2)*(f107A-71.5) + B1(3)*(f107-f107A+3.9)
     endif
     if (fexvir > 0.001) THEN
        r2 = fexvir
     else
        r2 =  B2(1) + B2(2)*(f107A-71.5) + B2(3)*(f107-f107A+3.9)
     endif
     do N = 13, Num_WaveLengths_High
        Solar_Flux(N) = (RFLUX(N) + ((r1-1.)*SCALE1(N)              &
             +  (r2-1.)*SCALE2(N)) / 1000.)
     enddo

  case (Hinteregger_Linear_Interp)

     ! do nothing

  case (Tobiska_EUV91)

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

     do N=16,55
        Solar_Flux(N) = TCHR0(N)        + &
             TCHR1(N)*hlymod + &
             TCHR2(N)*heimod + &
             TCOR0(N)        + &
             TCOR1(N)*f107   + &
             TCOR2(N)*f107A
     enddo

  case (WoodsAndRottman_10Nov88)

     DO N=15,55
        Solar_Flux(N) = WAR1(N)
     enddo

  case (WoodsAndRottman_20Jun89)

     DO N=15,55
        Solar_Flux(N) = WAR2(N)
     enddo

  end select

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

  do N=1,Num_WaveLengths_High

     IF (Solar_Flux(N) < 0.0) Solar_Flux(N) = 0.0

     IF ((WAVEL(N) < 251.0) .AND. (WAVES(N) > 15.0)) then
        Solar_Flux(N) = Solar_Flux(N)*xuvf
     endif

     !
     ! Convert to photons/m^2/s
     !

     Solar_Flux(N) = Solar_Flux(N) * 1.E9 * 10000.0

     !
     ! Convert to eV/m^2/s
     !


     wavelength_ave = (WAVEL(N) + WAVES(N))/2.0
     PhotonEnergy(N)= 6.626e-34*2.998e8/(wavelength_ave*1.0e-10)

     !
     !       Solar_Flux(N) = Solar_Flux(N) * &
     !            (Planck_Constant * Speed_Light) / &
     !            (wavelength_ave * 1.0e-10 * Element_Charge)

  enddo

  do N = 16,Num_WaveLengths_High
     Flux_of_EUV(N) = Solar_Flux(N)
  enddo

  do N = 1,15
     Flux_of_EUV(N) = Solar_Flux(N)
  enddo

  do N=1,Num_WaveLengths_Low
     NN = N+15
     Flux_of_EUV(NN) = 0.5*(EUV_Flux(N)+Solar_Flux(NN))
  enddo


  Flux_of_EUV = Flux_of_EUV/(SunOrbitEccentricity**2)


  do N=1,Num_WaveLengths_High
     wvavg(N)=(WAVEL(N)+WAVES(N))/2.
  enddo

 if (UseEUVData) then

    call start_timing("new_euv")
    SeeTime(:) = 0

    do N = 1, nSeeTimes
       SeeTime(N) = TimeSee(N)
    enddo

    tDiff = CurrentTime - SeeTime
    
    where (tDiff .lt. 0) tDiff = 1.e20
    iMin = minloc(tDiff)
    
    Timed_Flux = SeeFlux(:,iMin(1))

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
                Timed_Flux = SeeFlux(:,FlareStartIndex)
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
    
    !!need to convert from W/m^2 to photons/m^2/s
    do N=1,Num_WaveLengths_High 
        Flux_of_EUV(N) = Timed_Flux(N)*wvavg(N)*1.0e-10/(6.626e-34*2.998e8) &
             /(SunOrbitEccentricity**2)
     enddo
     call end_timing("new_euv")

  endif

  ! Second Spectra, provided by Steve Bougher....

!  do n = 1, nS2WaveLengths
!     S2PhotonEnergy(N) = 6.626e-34*2.998e8/(S2WaveLengths(n)*1.0e-10)
!  enddo

end subroutine calc_scaled_euv

 !-------------------------------------------------------------------
! Subroutine for initializing photoabsorption, photoionization,
! and branching ratio quantities.
!-------------------------------------------------------------------

subroutine init_euv

  use ModEUV

  implicit none

  integer :: N, NN

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
     PhotoAbs_CO2      = PhotoAbs_CO2        / 10000.0
     PhotoAbs_CO      = PhotoAbs_CO        / 10000.0

     PhotoIon_O2      = PhotoIon_O2        / 10000.0
     PhotoIon_CO2      = PhotoIon_CO2        / 10000.0
     PhotoIon_CO      = PhotoIon_CO        / 10000.0
     PhotoIon_OPlus4S = PhotoIon_OPlus4S   / 10000.0
     PhotoIon_N2      = PhotoIon_N2        / 10000.0
     PhotoIon_N       = PhotoIon_N         / 10000.0
     PhotoIon_OPlus2D = PhotoIon_OPlus2D   / 10000.0
     PhotoIon_OPlus2P = PhotoIon_OPlus2P   / 10000.0
     


!  do n = 1, nS2WaveLengths
!
!     S2PhotoAbsCO  = S2PhotoAbsCO * 1.0e-18 / 10000.0
!     S2PhotoAbsCO2 = S2PhotoAbsCO2 * 1.0E-18 / 10000.0
!     S2PhotoAbsN2  = S2PhotoAbsN2 * 1.0E-18 / 10000.0
!     S2PhotoAbsO2  = S2PhotoAbsO2 * 1.0E-18 / 10000.0
!     S2PhotoIonN2  = S2PhotoIonN2 * 1.0E-18 / 10000.0
!     S2PhotoIonCO  = S2PhotoIonCO * 1.0E-18 / 10000.0
!     S2PhotoIonCO2 = S2PhotoIonCO2 * 1.0E-18 / 10000.0
!     S2PhotoIonO2  = S2PhotoIonO2 * 1.0E-18 / 10000.0
!
!  enddo

end subroutine init_euv

subroutine Set_Euv(iError, StartTime, EndTime)

  use ModKind
  use ModEUV
  use ModInputs

  implicit none

  integer, intent(out)     :: iError
  real(Real8_), intent(in) :: EndTime, StartTime

  character(len=20)                       :: line, cline
  integer, dimension(7)                   :: TimeOfFlare,TimeArray
  real, dimension(6+Num_Wavelengths_High) :: temp

  real (Real8_) :: TimeDelay, BufferTime = 180.0

  logical :: NotDone
  integer :: i, iline, ioerror, nline, nILine = 1, iErrortemp = 0

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

end subroutine Set_Euv



  
  
  
