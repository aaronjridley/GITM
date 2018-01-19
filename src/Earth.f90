!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!-----------------------------------------------------------------------------
! $Id: Earth.f90,v 1.21 2017/08/05 01:33:09 ridley Exp $
!
! Author: Aaron Ridley, UMichigan
!
! Modified: AGB Oct 2013 - Corrected spelling of photoelectron heating
!                          efficiency variable
!-----------------------------------------------------------------------------

subroutine fill_photo(photoion, photoabs, photodis)

  use ModPlanet
  use ModEUV

  implicit none

  real, intent(out) :: photoion(Num_WaveLengths_High, nIons-1)
  real, intent(out) :: photoabs(Num_WaveLengths_High, nSpeciesTotal)
  real, intent(out) :: photodis(Num_WaveLengths_High, nSpeciesTotal)

  integer :: iSpecies, iWave

  PhotoAbs = 0.0
  PhotoIon = 0.0
  PhotoDis = 0.0

  photoabs(:,iO_3P_)  = PhotoAbs_O
  photoabs(:,iO2_)    = PhotoAbs_O2

  if (nSpecies > 2) then
     iSpecies = iN2_
     photoabs(:,iSpecies)    = PhotoAbs_N2
  endif
  if (nSpecies > 3) then
     iSpecies = iN_4S_
     photoabs(:,min(iSpecies,nSpecies))    = PhotoIon_N
  endif

  ! JMB:  06/25/2016
  if (nSpecies > 5) then
     iSpecies = iHe_
     photoabs(:,min(iSpecies,nSpecies))    = PhotoAbs_He
  endif


  ! This may need to be as defined below....
  photoion(:,iN2P_)   = PhotoIon_N2
  photoion(:,iO2P_)   = PhotoIon_O2
  photoion(:,iNP_)    = PhotoIon_N
  photoion(:,iO_4SP_) = PhotoIon_OPlus4S
  photoion(:,iO_2DP_) = PhotoIon_OPlus2D
  photoion(:,iO_2PP_) = PhotoIon_OPlus2P
  photoion(:,iHeP_)   = PhotoAbs_He

  do iWave = 1, Num_WaveLengths_High
     if (waves(iWave) >= 1250.0 .and. wavel(iWave) <= 1750.0) then
        PhotoDis(iWave, iO2_) = &
             photoabs(iWave,iO2_) - PhotoIon(iWave, iO2P_)
     endif

     if (waves(iWave) >= 800.0 .and. wavel(iWave) <= 1250.0) then
        PhotoDis(iWave, iN2_) = &
             photoabs(iWave,iN2_) - PhotoIon(iWave, iN2P_)
     endif

  enddo

  ! PE Ratio:  N2 + e- -> N(2D) + N(4S)
  pelecratio_N2(:,1) = PhotoElec_N2_N4S
  ! PE Ratio:  N2 + e- -> N+ + N(4S)
  pelecratio_N2(:,2) = PhotoElec_N2_NPlus
  ! PE Ratio:  N2 + e- -> N2+
  pelecratio_N2(:,3) = PhotoElec_N2_N2Plus

  ! PE Ratio:  O2 + e- -> O(4S) + O(3P)
  pelecratio_O2(:,1) = PhotoElec_O2_O3P
  ! PE Ratio:  O2 + e- -> O2+
  pelecratio_O2(:,2) = PhotoElec_O2_O2Plus
  ! PE Ratio:  O2 + e- -> O(4S)+ + O(3P)
  pelecratio_O2(:,3) = PhotoElec_O2_OPlus

end subroutine fill_photo

subroutine calc_planet_sources(iBlock)

  use ModInputs
  use ModSources
  use ModEUV
  use ModGITM
  use ModTime
  
  implicit none

  integer, intent(in) :: iBlock

  integer :: iAlt, iError, iDir, iLat, iLon

  real :: tmp2(nLons, nLats, nAlts)
  real :: tmp3(nLons, nLats, nAlts)
  real :: Omega(nLons, nLats, nAlts)
  real :: CO2Cooling(nLons, nLats, nAlts)

  LowAtmosRadRate = 0.0

  !\
  ! Cooling ----------------------------------------------------------
  !/

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 4) write(*,*) "=====> NO cooling", iproc, UseNOCooling

  call calc_co2(iBlock)

  if (UseCO2Cooling) then

     ! The 0.165 is derived from the TIEGCM (2.65e-13 / 1.602e-12)
     ! multiplied by 1e6 for /cm2 to /m2
     CO2Cooling = 0.165e6 * NDensityS(1:nLons,1:nLats,1:nAlts,iCO2_,iBlock)*&
          exp(-960.0/( &
            Temperature(1:nLons,1:nLats,1:nAlts,iBlock)* &
            TempUnit(1:nLons,1:nLats,1:nAlts))) * &
          MeanMajorMass(1:nLons,1:nLats,1:nAlts) * ( &
           (NDensityS(1:nLons,1:nLats,1:nAlts,iO2_,iBlock)/Mass(iO2_) + &
            NDensityS(1:nLons,1:nLats,1:nAlts,iN2_,iBlock)/Mass(iN2_)) * &
           2.5e-15 / 1e6 + &
           (NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock)/Mass(iO_3P_)) * &
           1.0e-12 / 1e6) * 1.602e-19

  else
     CO2Cooling = 0.0
  endif

  if (UseNOCooling) then

     !  [NO] cooling 
     ! [Reference: Kockarts,G., G.R.L.,VOL.7, PP.137-140,Feberary 1980 ]
 
     Omega = 3.6e-17 * NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock) /      &
          (3.6e-17 * NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock) + 13.3)

     ! We need to check this out. I don't like the first / sign....

     NOCooling = Planck_Constant * Speed_Light / &
          5.3e-6 * &
          Omega * 13.3 *  &
          exp(- Planck_Constant * Speed_Light / &
          (5.3e-6 * Boltzmanns_Constant * &
          Temperature(1:nLons,1:nLats,1:nAlts,iBlock)* &
          TempUnit(1:nLons,1:nLats,1:nAlts))) * &
          NDensityS(1:nLons,1:nLats,1:nAlts,iNO_,iBlock)

     RadiativeCooling2d = 0.0
     do iAlt=1,nAlts
        RadiativeCooling2d(1:nLons, 1:nLats) = &
             RadiativeCooling2d(1:nLons, 1:nLats) + &
             NOCooling(1:nLons,1:nLats,iAlt) * dAlt_GB(1:nLons,1:nLats,iAlt,iBlock)
     enddo

     NOCooling = NOCooling / TempUnit(1:nLons,1:nLats,1:nAlts) / &
          (Rho(1:nLons,1:nLats,1:nAlts,iBlock)*cp(:,:,1:nAlts,iBlock))

  else

     NOCooling = 0.0

  endif

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 4) write(*,*) "=====> UseOCooling", iproc, UseOCooling

  if (UseOCooling) then 

     ! [O] cooling 
     ! Reference: Kockarts, G., Plant. Space Sci., Vol. 18, pp. 271-285, 1970
     ! We reduce the LTE 63-um cooling rate by a factor of 2 for 
     ! the non-LTE effects.[Roble,1987]         

     tmp2 = exp(-228./(Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*&
          TempUnit(1:nLons,1:nLats,1:nAlts)))
     tmp3 = exp(-326./(Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*&
          TempUnit(1:nLons,1:nLats,1:nAlts)))

     ! In erg/cm3/s
     OCooling = (1.69e-18*tmp2 + 4.59e-20*tmp3) * &
          (NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock)/1.0e6) / &
          (1.0 + 0.6*tmp2 + 0.2*tmp3)
     ! In w/m3/3
     OCooling = OCooling/10.0
     ! In our special units:
     OCooling = OCooling/ TempUnit(1:nLons,1:nLats,1:nAlts) / &
          (Rho(1:nLons,1:nLats,1:nAlts,iBlock)*cp(:,:,1:nAlts,iBlock))

  else

     OCooling = 0.0

  endif

!  do iAlt = 1,15
!     write(*,*) 'no, co2 : ',iAlt, Altitude_GB(1,1,iAlt,iBlock)/1e3, &
!          NOCooling(1,1,iAlt), CO2Cooling(1,1,iAlt)
!  enddo

  RadCooling(1:nLons,1:nLats,1:nAlts,iBlock) = &
       OCooling + NOCooling + CO2Cooling


  PhotoElectronHeating(:,:,:,iBlock) = 0.0
  PhotoElectronHeating(:,:,:,iBlock) = &
       PhotoElectronHeatingEfficiency * &
       35.0*1.602e-19*&
       ( &
       EuvIonRateS(:,:,:,iO2P_,iBlock)* &
       NDensityS(1:nLons,1:nLats,1:nAlts,iO2_,iBlock) + &
       EuvIonRateS(:,:,:,iN2P_,iBlock)* &
       NDensityS(1:nLons,1:nLats,1:nAlts,iN2_,iBlock) + &
       EuvIonRateS(:,:,:,iO_4SP_,iBlock)* &
       NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock) + &
       EuvIonRateS(:,:,:,iO_2DP_,iBlock)* &
       NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock) + &
       EuvIonRateS(:,:,:,iO_2PP_,iBlock)* &
       NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock))
  
  PhotoElectronHeating(:,:,:,iBlock) = &
       PhotoElectronHeating(:,:,:,iBlock) / &
       Rho(1:nLons,1:nLats,1:nAlts,iBlock) / &
       cp(1:nLons,1:nLats,1:nAlts,iBlock) / &
       TempUnit(1:nLons,1:nLats,1:nAlts)

!--------------------------------------------------------------------
! GLOW
!--------------------------------------------------------------------

if (UseGlow) then
     if (dt < 10000.) then
        if  (floor((tSimulation-dt)/DtGlow) /= &
             floor(tsimulation/DtGlow)) then   

           call start_timing("glow")
           isInitialGlow = .True.

           if (iDebugLevel > 4) write(*,*) "=====> going into get_glow", iproc

           do iLat = 1, nLats
              do iLon = 1, nLons

                 call get_glow(iLon,iLat,iBlock)
                 
              enddo
           enddo

           call end_timing("glow")

        endif
     endif
     PhotoElectronDensity(:,:,:,:,iBlock) = PhotoElectronRate(:,:,:,:,iBlock) * dt
  endif


end subroutine calc_planet_sources

!---------------------------------------------------------------------
! Initialize Heating Rates
!---------------------------------------------------------------------

subroutine init_heating_efficiency

  use ModGITM, only: nLons, nLats, nAlts, nBlocks, Altitude_GB
  use ModEUV, only: HeatingEfficiency_CB, eHeatingEfficiency_CB
  use ModInputs, only: NeutralHeatingEfficiency

  implicit none

  integer :: iLon, iLat, iAlt
  !------------------------------------------------------------------
  HeatingEfficiency_CB(:,:,:,1:nBlocks) = NeutralHeatingEfficiency
!  max(0.1, &
!       0.40 - &
!       5.56e-5*(Altitude_GB(1:nLons,1:nLats,1:nAlts,1:nBlocks)/1000 - 165)**2)

  where(Altitude_GB(1:nLons,1:nLats,1:nAlts,1:nBlocks)/1000. > 150.)
     eHeatingEfficiency_CB(:,:,:,1:nBlocks) = 0.04
!!! min(0.4, &
!!!     0.04 + &
!!!     0.05*(Altitude_GB(1:nLons,1:nLats,1:nAlts,1:nBlocks)/1000 - 150)/100)
  elsewhere        
     eHeatingEfficiency_CB(:,:,:,1:nBlocks) = max(0.000001, &
          0.05 + &
          0.07*(Altitude_GB(1:nLons,1:nLats,1:nAlts,1:nBlocks)/1000 - 200)/100)
  end where

end subroutine init_heating_efficiency

!---------------------------------------------------------------------
! Calculate Eddy Diffusion Coefficient
!---------------------------------------------------------------------

subroutine calc_eddy_diffusion_coefficient(iBlock)

  use ModSizeGITM
  use ModGITM, only: pressure
  use ModInputs, only: EddyDiffusionPressure0,EddyDiffusionPressure1, &
       EddyDiffusionCoef
  use ModSources, only: KappaEddyDiffusion

  implicit none

  integer, intent(in) :: iBlock
  integer :: iAlt, iLat, iLon

  KappaEddyDiffusion=0.
  do iAlt = -1, nAlts+2

     do iLat = 1, nLats
        do iLon = 1, nLons

           if (pressure(iLon,iLat,iAlt,iBlock) >EddyDiffusionPressure0) then
              KappaEddyDiffusion(iLon,iLat,iAlt,iBlock) = EddyDiffusionCoef
              
           else if (pressure(iLon,iLat,iAlt,iBlock) > &
                EddyDiffusionPressure1) then

              KappaEddyDiffusion(iLon,iLat,iAlt,iBlock) = EddyDiffusionCoef * &
                   (pressure(iLon,iLat,iAlt,iBlock) - &
                   EddyDiffusionPressure1)/&
                   (EddyDiffusionPressure0 - EddyDiffusionPressure1)

           endif
        enddo
     enddo
  enddo

end subroutine calc_eddy_diffusion_coefficient

subroutine set_planet_defaults

  use ModPlanet
  use ModInputs

  implicit none
  
  iNeutralDensityOutputList(iN_4S_)=.false.
  iNeutralDensityOutputList(iHe_)=.false.
  iNeutralDensityOutputList(iN_2D_)=.false.
  iNeutralDensityOutputList(iN_2P_)=.false.
  iNeutralDensityOutputList(iH_)=.false.
  iNeutralDensityOutputList(iCO2_)=.false.
  iNeutralDensityOutputList(iO_1D_)=.false.

  iIonDensityOutputList(iNP_)=.false.
  iIonDensityOutputList(iO_2DP_)=.false.
  iIonDensityOutputList(iO_2PP_)=.false.
  iIonDensityOutputList(iHP_)=.false.
  iIonDensityOutputList(iHeP_)=.false.

  iTemperatureOutputList(2)=.false.
  iTemperatureOutputList(3)=.false.
  
  return

end subroutine set_planet_defaults

