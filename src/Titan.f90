! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

!-----------------------------------------------------------------------------
! $Id: Titan.f90,v 1.3 2017/07/12 18:51:35 jmbell Exp $
!
! Author: Aaron Ridley, UMichigan
!
! Modified: AGB Oct 2013 - Corrected spelling of photoelectron heating
!                          efficiency variable
!-----------------------------------------------------------------------------

subroutine fill_photo(photoion, photoabs, photodis)

  use ModPlanet
  use ModInputs
  use ModEUV

  implicit none

  real, intent(out) :: photoion(Num_WaveLengths_High, nIons-1)
  real, intent(out) :: photoabs(Num_WaveLengths_High, nSpeciesTotal)
  real, intent(out) :: photodis(Num_WaveLengths_High, nSpeciesTotal)

  integer :: iSpecies, iWave, NWH, i, iIon
  NWH = Num_WaveLengths_High

  PhotoAbs = 0.0
  PhotoIon = 0.0
  PhotoDis = 0.0

  ! Net Absorption:  Used in Solar EUV Heating
  photoabs(1:NWH,iN2_)   = PhotoAbs_N2(1:NWH)
  photoabs(1:NWH,iCH4_)  = PhotoAbs_CH4(1:NWH)
  !photoabs(1:NWH,iHCN_)  = PhotoAbs_HCN(1:NWH)
  !photoabs(1:NWH,iH2_)   = PhotoAbs_H2(1:NWH)
  !photoabs(1:NWH,iC2H4_) = PhotoAbs_C2H4(1:NWH)

!  photoabs(1:NWH,iC2H2_) = PhotoAbs_C2H2(1:NWH)
!! Net N2 and CH4 Photo-dissociation 
   photodis(1:NWH,iN2_) = QuantumYield_N2_N4S(1:NWH) * &
                          PhotoAbs_N2(1:NWH)

   photodis(1:NWH,iCH4_)=  (   QuantumYield_CH4_CH3(1:NWH)  + &
                               QuantumYield_CH4_3CH2(1:NWH) + &
                               QuantumYield_CH4_1CH2(1:NWH) + &
                               QuantumYield_CH4_CH(1:NWH)  ) *  &
                               PhotoAbs_CH4(1:NWH)

  photodis(1:NWH,iH2_)= PhotoAbs_H2(1:NWH)
  photodis(1:NWH,iHCN_)= PhotoAbs_HCN(1:NWH)
  photodis(1:NWH,iC2H4_)= PhotoAbs_C2H4(1:NWH)

!--------
!   N2 Photolytic Reactions
  ! Rxn1:  hv + N2 -> N(2D) (65%) + N(4S)(35%)
  photodiss_N2(:,1)  = PhotoAbs_N2*QuantumYield_N2_N4S
  ! Rxn3:  hv + N2 -> N+ + N(4S)
  photodiss_N2(:,2)  = PhotoAbs_N2(:)*QuantumYield_N2_NPlus(:)
  ! Rxn3:  hv + N2 -> N2+ 
  photodiss_N2(:,3)  = PhotoAbs_N2(:)*QuantumYield_N2_N2Plus(:)

  ! Photoelectron contributions (Solomon and Qian [2005] for N2)
  ! PE Ratio:  N2 + e- -> N(2D) + N(4S)
  pelecratio_N2(:,1) = PhotoElec_N2_N4S
  ! PE Ratio:  N2 + e- -> N+ + N(4S)
  pelecratio_N2(:,2) = PhotoElec_N2_NPlus
  ! PE Ratio:  N2 + e- -> N2+
  pelecratio_N2(:,3) = PhotoElec_N2_N2Plus

!--------
! Photolysis of CH4 
!
!--------
!   CH4 Photolytic Reactions
  ! Rxn1:  hv + CH4 -> CH3 + H
  photodiss_CH4(:,1)  = PhotoAbs_CH4*QuantumYield_CH4_CH3
  ! Rxn3:  hv + CH4 -> 1CH2 + H2
  photodiss_CH4(:,2)  = PhotoAbs_CH4(:)*QuantumYield_CH4_1CH2(:)
  ! Rxn3:  hv + CH4 -> 3CH2 + H
  photodiss_CH4(:,3)  = PhotoAbs_CH4(:)*QuantumYield_CH4_3CH2(:)
  ! Rxn3:  hv + CH4 -> CH + H2 + H
  photodiss_CH4(:,4)  = PhotoAbs_CH4(:)*QuantumYield_CH4_CH(:)
  ! Rxn3:  hv + CH4 -> CH4+ 
  photodiss_CH4(:,5)  = PhotoAbs_CH4(:)*QuantumYield_CH4_CH4Plus(:)
  ! Rxn3:  hv + CH4 -> CH3+ + H 
  photodiss_CH4(:,6)  = PhotoAbs_CH4(:)*QuantumYield_CH4_CH3Plus(:)
  ! Rxn3:  hv + CH4 -> CH2+ + H 
  photodiss_CH4(:,7)  = PhotoAbs_CH4(:)*QuantumYield_CH4_CH2Plus(:)
  ! Rxn3:  hv + CH4 -> CH+ + H2 + H
  photodiss_CH4(:,8)  = PhotoAbs_CH4(:)*QuantumYield_CH4_CHPlus(:)
  ! Rxn3:  hv + CH4 -> C+ + 2*H2 
  photodiss_CH4(:,9)  = PhotoAbs_CH4(:)*QuantumYield_CH4_CPlus(:)
  ! Rxn3:  hv + CH4 -> 1CH2 + H2Plus
  photodiss_CH4(:,10)  = PhotoAbs_CH4(:)*QuantumYield_CH4_H2Plus(:)
  ! Rxn3:  hv + CH4 -> CH3 + HPlus
  photodiss_CH4(:,11)  = PhotoAbs_CH4(:)*QuantumYield_CH4_HPlus(:)

  ! Photoelectron contributions (Solomon and Qian [2005] for N2)
  ! PE Ratio:  CH4 + e* -> No Reactions
  pelecratio_CH4(:,1:4) = 0.0
  ! PE Ratio:  CH4 + e* -> CH4+ + e-
  pelecratio_CH4(:,5) = 5.0*PhotoElec_CH4_CH4Plus
  ! PE Ratio:  CH4 + e- -> CH3+ + H + e-
  pelecratio_CH4(:,6) = 20.0*PhotoElec_CH4_CH3Plus
  ! Ignore other ions
  pelecratio_CH4(:,7:11) = 0.0

! Old Method
  photodis(1:NWH,iCH3_)=  QuantumYield_CH4_CH3(1:NWH) * &
                          PhotoAbs_CH4(1:NWH)

  photodis(1:NWH,i3CH2_)= QuantumYield_CH4_3CH2(1:NWH) * &
                          PhotoAbs_CH4(1:NWH)

  photodis(1:NWH,i1CH2_)= QuantumYield_CH4_1CH2(1:NWH) * &
                          PhotoAbs_CH4(1:NWH)

  photodis(1:NWH,iCH_)= QuantumYield_CH4_CH(1:NWH) * &
                          PhotoAbs_CH4(1:NWH)

  photoion(1:NWH,iCH3P_)=   QuantumYield_CH4_CH3Plus(1:NWH)*  &
                            PhotoAbs_CH4(1:NWH)   

!!! C2H2
  ! Rxn15: hv + C2H2 -> C2H + H
  photodiss_C2H2(1:NWH,1) =  QuantumYield_C2H2_C2H_H(1:NWH)*  &
                            PhotoAbs_C2H2(1:NWH)   

  ! Rxn16: hv + C2H2 -> C2 + H2
  photodiss_C2H2(1:NWH,2) =  QuantumYield_C2H2_C2_H2(1:NWH)*  &
                            PhotoAbs_C2H2(1:NWH)   

  ! Rxn17: hv + C2H2 -> C2H2+ 
  photodiss_C2H2(1:NWH,3) =  QuantumYield_C2H2_C2H2Plus(1:NWH)*  &
                            PhotoAbs_C2H2(1:NWH)   

  ! Rxn18: hv + C2H2 -> C2H+ + H 
  photodiss_C2H2(1:NWH,4) =  QuantumYield_C2H2_C2HPlus_H(1:NWH)*  &
                            PhotoAbs_C2H2(1:NWH)   

  ! Rxn19: hv + C2H2 -> CH+ + CH 
  photodiss_C2H2(1:NWH,5) =  QuantumYield_C2H2_CHPlus_CH(1:NWH)*  &
                            PhotoAbs_C2H2(1:NWH)   

  ! Rxn19: hv + C2H2 -> C+ + CH2 
  photodiss_C2H2(1:NWH,6) =  QuantumYield_C2H2_CHPlus_CH(1:NWH)*  &
                            PhotoAbs_C2H2(1:NWH)   

  ! Rxn20: hv + C2H2 -> H+ + C2H
  photodiss_C2H2(1:NWH,7) =  QuantumYield_C2H2_HPlus_C2H(1:NWH)*  &
                            PhotoAbs_C2H2(1:NWH)   


!!! C2H4
  ! Rxn22: hv + C2H4 -> C2H2 + H2
  photodiss_C2H4(1:NWH,1) =  QuantumYield_C2H4_C2H2_H2(1:NWH)*  &
                            PhotoAbs_C2H4(1:NWH)   

  ! Rxn23: hv + C2H4 -> C2H2 + 2H
  photodiss_C2H4(1:NWH,2) =  QuantumYield_C2H4_C2H2_2H(1:NWH)*  &
                            PhotoAbs_C2H4(1:NWH)   

  ! Rxn24: hv + C2H4 -> C2H4+ 
  photodiss_C2H4(1:NWH,3) =  QuantumYield_C2H4_C2H4Plus(1:NWH)*  &
                            PhotoAbs_C2H4(1:NWH)   

  ! Rxn25: hv + C2H4 -> C2H3+ + H 
  photodiss_C2H4(1:NWH,4) =  QuantumYield_C2H4_C2H3Plus_H(1:NWH)*  &
                            PhotoAbs_C2H4(1:NWH)   

  ! Rxn26: hv + C2H4 -> C2H2+ + H2 
  photodiss_C2H4(1:NWH,5) =  QuantumYield_C2H4_C2H2Plus_H2(1:NWH)*  &
                            PhotoAbs_C2H4(1:NWH)   

  ! Rxn27: hv + C2H4 -> C2H2+ + H2 
  photodiss_C2H4(1:NWH,6) =  QuantumYield_C2H4_C2HPlus_H2_H(1:NWH)*  &
                            PhotoAbs_C2H4(1:NWH)   


!!! H Absorption
  ! Rxn28: hv + H -> H+
  photoion_H(1:NWH) = PhotoAbs_H(1:NWH)   

!!! H2 Absorption
!
!  ! Rxn29: H2 -> 2 H
!  photodis(1:NWH,iH_)= PhotoAbs_H2(1:NWH)*QuantumYield_H2_2H(1:NWH)
!  ! Rxn30: H2 ->  H2+
!  photoion(1:NWH,iH2P_)= PhotoAbs_H2(1:NWH)*QuantumYield_H2_H2P(1:NWH)
!  ! Rxn31: H2 ->  H + H+
!  photoion(1:NWH,iHP_)= PhotoAbs_H2(1:NWH)*QuantumYield_H2_H_HP(1:NWH)
!
  ! Rxn32: N + hv -> N+
  Newphotoion_N(1:NWH) = PhotoAbs_N(1:NWH)   
!!! HCN
  ! Rxn33: HCN + hv -> H + CN
  photodiss_HCN(1:NWH,1) =  QuantumYield_HCN_H_CN(1:NWH)*  &
                            PhotoAbs_HCN(1:NWH)   

  ! Rxn34: HCN + hv -> HCN+
  photodiss_HCN(1:NWH,2) =  QuantumYield_HCN_HCNP(1:NWH)*  &
                            PhotoAbs_HCN(1:NWH)   

end subroutine fill_photo

subroutine calc_planet_sources(iBlock)

  use ModInputs
  use ModSources
  use ModEUV
  use ModGITM
  use ModTime
  use ModPlanet, only : HCNCoolRatio
  
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

!
! Do Nothing
  if( (floor((tSimulation - Dt)/dTCooling) .eq. floor(tSimulation/dTCooling)) ) then
     if( (floor((tSimulation - Dt)/dTRatio) .eq. floor(tSimulation/dTRatio)) ) then
       ! do Nothing
     else
        !write(*,*) 'Call cool to space'
        call calc_cooltospace(iBlock)
     endif
  elseif(maxval(HCNCoolRatio) .eq. 0.0) then
     call calc_radcooling(iBlock)
  else
     call calc_radcooling(iBlock)
  endif
  

  RadCooling(1:nLons,1:nLats,1:nAlts,iBlock) = &
    RadCoolingRate(1:nLons,1:nLats,1:nAlts,iBlock)/&
      ( cp(1:nLons,1:nLats,1:nAlts,iBlock)*&
       rho(1:nLons,1:nLats,1:nAlts,iBlock)  )

  RadCooling(1:nLons,1:nLats,1:nAlts,iBlock) = &
     RadCooling(1:nLons,1:nLats,1:nAlts,iBlock)/&
       TempUnit(1:nLons,1:nLats,1:nAlts) 

! call calc_saturn_tides(iBlock)
  TideADep(:,:,:,:,iBlock) = 0.0
  TideAIndp(:,:,:,:,iBlock)= 0.0
!  RadCooling(1:nLons,1:nLats,1:nAlts,iBlock) = &
!    0.001

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 4) write(*,*) "=====> HCN cooling", iproc, UseNOCooling


  ! Residual Earth Stuff
!  call calc_co2(iBlock)
!
!  if (UseCO2Cooling) then
!
!     ! The 0.165 is derived from the TIEGCM (2.65e-13 / 1.602e-12)
!     ! multiplied by 1e6 for /cm2 to /m2
!     CO2Cooling = 0.165e6 * NDensityS(1:nLons,1:nLats,1:nAlts,iCO2_,iBlock)*&
!          exp(-960.0/( &
!            Temperature(1:nLons,1:nLats,1:nAlts,iBlock)* &
!            TempUnit(1:nLons,1:nLats,1:nAlts))) * &
!          MeanMajorMass(1:nLons,1:nLats,1:nAlts) * ( &
!           (NDensityS(1:nLons,1:nLats,1:nAlts,iO2_,iBlock)/Mass(iO2_) + &
!            NDensityS(1:nLons,1:nLats,1:nAlts,iN2_,iBlock)/Mass(iN2_)) * &
!           2.5e-15 / 1e6 + &
!           (NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock)/Mass(iO_3P_)) * &
!           1.0e-12 / 1e6) * 1.602e-19
!
!  else
!     CO2Cooling = 0.0
!  endif
!
!  if (UseNOCooling) then
!
!     !  [NO] cooling 
!     ! [Reference: Kockarts,G., G.R.L.,VOL.7, PP.137-140,Feberary 1980 ]
! 
!     Omega = 3.6e-17 * NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock) /      &
!          (3.6e-17 * NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock) + 13.3)
!
!     ! We need to check this out. I don't like the first / sign....
!
!     NOCooling = Planck_Constant * Speed_Light / &
!          5.3e-6 * &
!          Omega * 13.3 *  &
!          exp(- Planck_Constant * Speed_Light / &
!          (5.3e-6 * Boltzmanns_Constant * &
!          Temperature(1:nLons,1:nLats,1:nAlts,iBlock)* &
!          TempUnit(1:nLons,1:nLats,1:nAlts))) * &
!          NDensityS(1:nLons,1:nLats,1:nAlts,iNO_,iBlock)
!
!     NOCooling = NOCooling / TempUnit(1:nLons,1:nLats,1:nAlts) / &
!          (Rho(1:nLons,1:nLats,1:nAlts,iBlock)*cp(:,:,1:nAlts,iBlock))
!
!  else
!
!     NOCooling = 0.0
!
!  endif
!
!  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
!  if (iDebugLevel > 4) write(*,*) "=====> UseOCooling", iproc, UseOCooling
!
!  if (UseOCooling) then 
!
!     ! [O] cooling 
!     ! Reference: Kockarts, G., Plant. Space Sci., Vol. 18, pp. 271-285, 1970
!     ! We reduce the LTE 63-um cooling rate by a factor of 2 for 
!     ! the non-LTE effects.[Roble,1987]         
!
!     tmp2 = exp(-228./(Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*&
!          TempUnit(1:nLons,1:nLats,1:nAlts)))
!     tmp3 = exp(-326./(Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*&
!          TempUnit(1:nLons,1:nLats,1:nAlts)))
!
!     ! In erg/cm3/s
!     OCooling = (1.69e-18*tmp2 + 4.59e-20*tmp3) * &
!          (NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock)/1.0e6) / &
!          (1.0 + 0.6*tmp2 + 0.2*tmp3)
!     ! In w/m3/3
!     OCooling = OCooling/10.0
!     ! In our special units:
!     OCooling = OCooling/ TempUnit(1:nLons,1:nLats,1:nAlts) / &
!          (Rho(1:nLons,1:nLats,1:nAlts,iBlock)*cp(:,:,1:nAlts,iBlock))
!
!  else
!
!     OCooling = 0.0
!
!  endif

!  do iAlt = 1,15
!     write(*,*) 'no, co2 : ',iAlt, Altitude_GB(1,1,iAlt,iBlock)/1e3, &
!          NOCooling(1,1,iAlt), CO2Cooling(1,1,iAlt)
!  enddo

!  RadCooling(1:nLons,1:nLats,1:nAlts,iBlock) = &
!       OCooling + NOCooling + CO2Cooling


! Note: Photoelectron heating is basically
!       Swartz and Nisbet [1972].
! This could be updated with Smithro and Solomon
!  PhotoElectronHeating(:,:,:,iBlock) = 0.0
!  PhotoElectronHeating(:,:,:,iBlock) = &
!       PhotoElectronHeatingEfficiency * &
!       35.0*1.602e-19*&
!       ( &
!       EuvIonRateS(:,:,:,iO2P_,iBlock)* &
!       NDensityS(1:nLons,1:nLats,1:nAlts,iO2_,iBlock) + &
!       EuvIonRateS(:,:,:,iN2P_,iBlock)* &
!       NDensityS(1:nLons,1:nLats,1:nAlts,iN2_,iBlock) + &
!       EuvIonRateS(:,:,:,iO_4SP_,iBlock)* &
!       NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock) + &
!       EuvIonRateS(:,:,:,iO_2DP_,iBlock)* &
!       NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock) + &
!       EuvIonRateS(:,:,:,iO_2PP_,iBlock)* &
!       NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock))
!  
!  PhotoElectronHeating(:,:,:,iBlock) = &
!       PhotoElectronHeating(:,:,:,iBlock) / &
!       Rho(1:nLons,1:nLats,1:nAlts,iBlock) / &
!       cp(1:nLons,1:nLats,1:nAlts,iBlock) / &
!       TempUnit(1:nLons,1:nLats,1:nAlts)

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
  HeatingEfficiency_CB(:,:,:,1:nBlocks) = 0.22

end subroutine init_heating_efficiency

!---------------------------------------------------------------------
! Calculate Eddy Diffusion Coefficient
!---------------------------------------------------------------------

subroutine calc_eddy_diffusion_coefficient(iBlock)

  use ModSizeGITM
  use ModGITM, only: NDensity, Pressure
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

           KappaEddyDiffusion(iLon,iLat,iAlt,iBlock) = &
            EddyDiffusionCoef*&
            sqrt(NDensity(iLon,iLat,1,iBlock)/NDensity(iLon,iLat,iAlt,iBlock))

!           if (pressure(iLon,iLat,iAlt,iBlock) >EddyDiffusionPressure0) then
!              KappaEddyDiffusion(iLon,iLat,iAlt,iBlock) = EddyDiffusionCoef
!              
!           else if (pressure(iLon,iLat,iAlt,iBlock) > &
!                EddyDiffusionPressure1) then
!
!              KappaEddyDiffusion(iLon,iLat,iAlt,iBlock) = EddyDiffusionCoef * &
!                   (pressure(iLon,iLat,iAlt,iBlock) - &
!                   EddyDiffusionPressure1)/&
!                   (EddyDiffusionPressure0 - EddyDiffusionPressure1)
!
!           endif

        enddo
     enddo
  enddo

end subroutine calc_eddy_diffusion_coefficient

subroutine set_planet_defaults

  use ModInputs
end subroutine !set_planet_defaults


      subroutine calc_radcooling(iBlock)
!  ---------------------------------------------------------------------
!  Purpose: 
!  =======
!  This is the Cooling Routine to be used for the Titan TGCM
!  This version incorporates all heating and cooling terms into
!  the full line by line radiative transfer calculations.
!     Date            Programmer                   Version
!  ==========     ==================            ==============
!  09-08-04           Jared Bell                 Original(1.0) 
!  12-02-04           Jared Bell                 Update (1.1)   
!  12-07-04           Jared Bell                 Update (1.2)   
!  01-27-05           Jared Bell                 Consistency Check   
!  01-27-05           Jared Bell                 Added Heating 1 
!  01-27-05           Jared Bell                 Added Heating 2 
!  01-28-05           Jared Bell                 Added Dynamic Memory 
!  09-07-05           Jared Bell                 Modify to use TAU 
!  09-12-05           Jared Bell                 Utilized Linear Interp. 
!  09-16-05           Jared Bell                 Expint more accurate 
!                                                and more efficient 
!  03-02-06           Jared Bell                 GITM                 
!  0l-07-07           Jared Bell                 Computational Efficiency
!                                                Update 
!  05-28-11           Jared Bell                 Removed all 5-D Arrays 
!                                                (nLon, nLat, nAlts, nLines, nFreqas)
!                                                Implemented Loops Instead
!  -----------------------------------------------------------------

      use ModPlanet
      use ModGITM
      use ModConstants, only:  Boltzmanns_Constant, Speed_Light, &
                               Planck_Constant
      use ModInputs, only:  iDebugLevel
      use ModSources, only: RadCoolingRate
      use ModTime, only: tSimulation, CurrentTime
      use ModIndicesInterfaces

      implicit none

! INPUTS:----------------------------------------------------------

      integer,intent(in):: iBlock 

! LOCAL VARS:------------------------------------------------------

      integer :: i,j,l,m,k,n,ii,tj 
      integer :: iLat, iLon, iAlt, iLine, iFreq

      integer, parameter :: NLVLS = nAlts+4

      integer, parameter :: NLINES = rotlines 

!      integer, parameter :: NFREQ = rotfreqs 
      integer, parameter :: NewNFREQ = newrotfreqs 
      integer, parameter :: NPTS = rotpts 

!      integer, parameter :: blon = 1 
!      integer, parameter :: elon = nLons 
!      integer, parameter :: blat = 1 
!      integer, parameter :: elat = nLats 
!
!      real,dimension(-1:nAlts+2) :: &
!       TCOOLa

      real ::   &
       FTH

      real,dimension(-1:nAlts+2) ::    &
       T,       & 
       NHCN,    &
       fcp,     &
       COOL,    &
       HEAT1,   &
       HEAT4

      real,dimension(1:NLINES,-1:nAlts+2) ::  &
       gcool,   &
       gheat1,  &
       gheat4,  &
       INTEN,   &
       ALPHAD,  &
       ALPHAV,  &
       ALPHAL,  &    ! VOIGT
       DELTAV,  &
       VULIM,   &
       VLLIM,   &
       VJM,     &
       TM,      &
       ZM,      &
       VTH,     &
       NHCNM,   &
       GammaNew,  &   ! VOIGT
       YRATIO         ! VOIGT

!!!! FWHM of the Doppler and Lorentz Profiles
      real,dimension(1:NLINES,-1:nAlts+2) ::  &
        FWHML, FWHMD, FWHMV

!!! Coefficients for Combining Gaussian and Lorentz shapes
!!! DRatio used by Liu et al. [2001] in JQSRT
      real :: DRatio, CoefL, CoefG
!      real,dimension(1:NLINES,-1:nAlts+2) ::  &
!        CoefL, CoefG

      real,dimension(1:NLINES,-1:nAlts+2) ::  &
       NEW, &
       OLD, &
       Y1,  &
       Y2,  &
       Y

      real,dimension(1:NLINES,-1:nAlts+2) ::  &
       PRODCOOL,   &
       PRODHEAT1,  &
       PRODHEAT4

      real,dimension(-1:nAlts+2) ::    &
       SMCOOL,   &
       SMHEAT1,  &
       SMHEAT4
!
!

!      real,dimension(-1:nAlts+2) ::  &
!       NHCNT 
!   
! --------------BEGIN TAU VARS--------------------------------------- 
! The Following Variables are used solely within the Tau Integration
! 
      real :: DFACT
      real,dimension(1:NewNFREQ,1:NLINES,-1:nAlts+2) ::  &
       NewTAU,  &
       SumTAU,  &
       NewDOP,  &  ! the DOPPLER Lineshape
       NewLorentz,  &  ! the DOPPLER Lineshape
       NewX,  &
       NewV, &
       NewYV,  &   ! VOIGT
       NewVOI,  &  ! the VOIGT Lineshape
       NewDopFactor

      real,dimension(1:NewNFREQ,1:NLINES,-1:nAlts+2) ::  &
       NewBv,  &
       NewBalpha

      real :: &
       Bfactor

      real :: timestart, timeend

      real :: &
       NONLTEFACTOR, &
       NONLTEPHI, &
       NONLTEQUENCHING, &
       NONLTEEINSTEIN

      integer :: MiddleFreqIndex, FreqIndexMax

      real ::  HCNTopsideScaleHeight

      real    :: RadCoolWidth, RadCoolFlagValue
      real    :: DeltaAlt, ExpArg
      real    :: ScaleRadCool(1:nAlts)
      integer :: iAltFlag
      logical :: IsFound

      real    :: TestX, TestAns
      IsFound = .false.

   MiddleFreqIndex = 5
   FreqIndexMax = 4

  call start_timing("calc_radcooling")

      HEAT1 = 0.0D0
      HEAT4 = 0.0D0

     do iLon = 1, nLons
       do iLat = 1, nLats

       T(-1:nAlts+2) = &
        Temperature(iLon,iLat,-1:nAlts+2,iBlock)*TempUnit(iLon,iLat,-1:nAlts+2)

       NHCN(-1:nAlts+2) = &
        NDensityS(iLon,iLat,-1:nAlts+2,iHCN_,iBlock)*1.0e-06

       HCNTopsideScaleHeight = &
           -1.0*(T(nAlts+2)*Boltzmanns_Constant/&
              (Mass(iHCN_)*Gravity_GB(iLon,iLat,nAlts+2,iBlock)))


! ------------------------------------------------------------------------
! INTERPOLATING INTENSITY DATA TO TEMPERATURE GRID -----------------------
! ------------------------------------------------------------------------
             do iAlt = -1, nAlts+2
              do iLine = 1, NLINES
!\
! Speed_Light*100.0:   Converts cm^-1 --> Hz 
! 1.0e-04          : converts Sj from Hz*cm^2/molecule -> Hz*m^2/molecule
                INTEN(iLine,iAlt) =             &
                  ( Speed_Light*100.0)*(1.0e-4)*          &
                  exp( pmat(iLine,1)                         &
                  + pmat(iLine,2)*(T(iAlt)   )  &
                  + pmat(iLine,3)*(T(iAlt)**2)  &
                  + pmat(iLine,4)*(T(iAlt)**3)  &
                  + pmat(iLine,5)*(T(iAlt)**4)  &
                  + pmat(iLine,6)*(T(iAlt)**5)  &
                  + pmat(iLine,7)*(T(iAlt)**6)  &
                  + pmat(iLine,8)*(T(iAlt)**7)  & 
                  + pmat(iLine,9)*(T(iAlt)**8)  )    


              enddo! iLine = 1, NLINES
             enddo !iAlt = -1, nAlts+2

! ------------------------------------------------------------------------
! Frequency Integration Variables TM,VJM,ZM,NHCNM (1:nlines,1:nAlts)
! ------------------------------------------------------------------------
             do iAlt = -1, nAlts+2 
!              write(*,*) '=====> iAlt =', iAlt
              do iLine = 1, NLINES  

               VJM(iLine,iAlt) = freqhz(iLine)        ! Lines in Hz
                ZM(iLine,iAlt) = Altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0   ! m -> km
                TM(iLine,iAlt) = T(iAlt)       ! Temp in  K
            NHCNM(iLine,iAlt) = NDensityS(iLon,iLat,iAlt,iHCN_,iBlock)*1.0E-06  ! m^-3 -> cm^-3

! These provide the Alpha_L in cm^-1 (we call this GammaNew)
          GammaNew(iLine,iAlt) = &
                    ALPHAL0(iLine)*&
             (1.0 - NDensityS(iLon,iLat,iAlt,iHCN_,iBlock)/NDensity(iLon,iLat,iAlt,iBlock))*&
             (Pressure(iLon,iLat,iAlt,iBlock)/(101325.00))*(298.0/T(iAlt))**ALPHALExp(iLine) + &
                    ALPHAS0(iLine)*&
             (NDensityS(iLon,iLat,iAlt,iHCN_,iBlock)/NDensity(iLon,iLat,iAlt,iBlock))*&
             (Pressure(iLon,iLat,iAlt,iBlock)/(101325.00))*(298.0/T(iAlt))**ALPHALExp(iLine)

!! AlphaL is the Half-width of the Lorentz profile (pressure-broadening) in Hz
          ALPHAL(iLine,iAlt) = &
                 GammaNew(iLine,iAlt)*(Speed_Light*100.0)   !! This is AlphaL in Hz

              enddo  ! End iLine
             enddo ! End iAlt

! ------------------------------------------------------------------------
! DOPPLER HALFWIDTHS AT GIVEN TEMPERATURES -------------------------------
! ------------------------------------------------------------------------
           FTH = SQRT(2.0*(Boltzmanns_Constant/Mass(iHCN_)))
           VTH = FTH*SQRT(TM)
           ALPHAD = (VJM*(VTH/(Speed_Light)))  !!! Our Doppler (Gaussian) half-width
           DELTAV = 5.0*ALPHAD
           VULIM = VJM + DELTAV
           VLLIM = VJM - DELTAV

! =================================================>
! Calculate the FWHM of the Convolved Voigt Profile
! Use the form taken from Olivier and Longbothom

           FWHMD = ALPHAD*(2.0*sqrt(2.0*alog(2.0)))  !! Our FUllWidth at Half-Max of Doppler
           FWHML = 2.0*ALPHAL

           FWHMV = 0.5346*FWHML + &
              sqrt(0.2166*FWHML**2.0 + FWHMD**2.0)

           ALPHAV = FWHMV/2.0  ! Our Voigt Half-width

!----------------
!WHILE LOOP BEGIN
!----------------
!
!      m = 1
!      do 
!      if (2**m .EQ. NFREQ) exit
!         m = m + 1
!      enddo
!----------------
!END WHILE LOOP 
!---------------

         NewV(MiddleFreqIndex,:,:) = &
                  (FreqQuadAbscissas(1)*(VULIM - VLLIM) + VULIM + VLLIM)*.5

       do k = 1 , FreqIndexMax
         NewV(MiddleFreqIndex - k,:,:)= &
              (FreqQuadAbscissas(k+1)*(VULIM - VLLIM) + VULIM + VLLIM)*.5
         NewV(MiddleFreqIndex + k,:,:)=&
              (-1.0*FreqQuadAbscissas(k+1)*(VULIM - VLLIM) + VULIM + VLLIM)*.5
       enddo

         DFACT = SQRT(alog(2.0)/PI)

        do iAlt = -1, nAlts+2
           do iLine = 1, NLINES
              ! Note that DRatio will vary 
              ! from -1.0 ( Pure Doppler)
              ! to 1.0 (Pure Lorentz)
              DRatio = (FWHML(iLine,iAlt) - FWHMD(iLine,iAlt))/&
                       (FWHML(iLine,iAlt) + FWHMD(iLine,iAlt))
 
              ! Coef L and CoefD determine the "relative mixing"
              ! Of our line shapes into the combined Voigt Shape
              CoefL = 0.68188 + 0.61293*DRatio - &
                      0.18384*DRatio**2.0 - &
                      0.11568*DRatio**3.0 

              CoefG = 0.32460 - 0.61825*DRatio + &
                      0.17681*DRatio**2.0 + &
                      0.12109*DRatio**3.0 

              do iFreq = 1, NewNFREQ

              ! Calculate the Doppler contribution
              
              NewX(iFreq,iLine,iAlt) =  &
                  SQRT( alog(2.0) )*  (  ( NewV(iFreq,iLine,iAlt)  - VJM(iLine,iAlt) )/ & 
                   ALPHAV(iLine,iAlt) ) 

              NewDOP(iFreq,iLine,iAlt) =  &
                    (DFACT/ALPHAV(iLine,iAlt))*    & 
                    EXP( -( NewX(iFreq,iLine,iAlt)**2.0) ) 

              NewYV(iFreq,iLine,iAlt) =  &
                    ( NewV(iFreq,iLine,iAlt)  - VJM(iLine,iAlt) )

              NewLorentz(iFreq,iLine,iAlt) =  &
                    (1.0/PI)*ALPHAV(iLine,iAlt)/    & 
                    ( NewYV(iFreq,iLine,iAlt)**2.0 + ALPHAV(iLine,iAlt)**2.0)

              NewVOI(iFreq,iLine,iAlt) = &
                   CoefL*NewLorentz(iFreq,iLine,iAlt) + &
                   CoefG*NewDOP(iFreq,iLine,iAlt) 

 !             write(*,*) '===============> iFreq', iFreq, NewVOI(iFreq,iLine,iAlt) 

              NewBalpha(iFreq,iLine,iAlt)  =             &
                    (  Planck_Constant/Boltzmanns_Constant)*  &
                    (  NewV(iFreq,iLine,iAlt) /        &
                            TM(iLine,iAlt)  )

              Bfactor = 1.0/ &
                        ( EXP(NewBalpha(iFreq,iLine,iAlt) ) - 1.0 ) 

              NewBv(iFreq,iLine,iAlt) =                   &
                     ( ( 2.0*Planck_Constant ) /        & 
                     ( (Speed_Light)**2 )  ) *                   &
                     ( NewV(iFreq,iLine,iAlt)**3 )*Bfactor

         enddo
       enddo
    enddo

! ------------------------------------------------------------------------
! Cooling CALCULATIONS
! Cool-To-Space Portion: 
! ------------------------------------------------------------------------

       OLD = 0.0
       NEW = 0.0
       NEW = NEW + FreqQuadWeights(1)*(NewBv(MiddleFreqIndex, :, :)*NewVOI(MiddleFreqIndex, :, :))

      do iFreq = 1 , FreqIndexMax
        Y1 = NewBv(MiddleFreqIndex - iFreq,:,:)*NewVOI(MiddleFreqIndex - iFreq,:,:)
        Y2 = NewBv(MiddleFreqIndex + iFreq,:,:)*NewVOI(MiddleFreqIndex + iFreq,:,:) 
        NEW = NEW + FreqQuadWeights(iFreq + 1)*(Y1 + Y2)
      enddo

      gcool(:,:) = (VULIM(:,:) - VLLIM(:,:) )*NEW*0.5
! Multiply by the Intensity of the lines and a geometric factor of 4pi
!
      PRODCOOL(1:NLINES,-1:nAlts+2) =  &
                       gcool(1:NLINES,-1:nAlts+2)* &
                       INTEN(1:NLINES,-1:nAlts+2)*(4.0D0*PI)
! Sum over the lines
      SMCOOL(-1:nAlts+2) = &
                SUM(PRODCOOL(1:NLINES,-1:nAlts+2),1)
! Multiply by HCN density 
      COOL(-1:nAlts+2)  =  &
                      SMCOOL(-1:nAlts+2)* &
                        NHCN(-1:nAlts+2)*1.0e+06   ! J/(m^3 s) 
! ------------------------------------------------------------------------
! END The Cool-To-Space Calculation: 
! ------------------------------------------------------------------------
      CALL Sumfgausstau(INTEN(1:NLINES,-1:nAlts+2), &
                        NHCNM(1:NLINES,-1:nAlts+2), &
                               ZM(1,-1:nAlts+2), &
               NewVOI(1:NewNFREQ,1:NLINES,-1:nAlts+2), &
               SumTAU(1:NewNFREQ,1:NLINES,-1:nAlts+2), &
              HCNTopsideScaleHeight )

      NewTAU = SumTAU

!
       CALL Newheat1gauss(NewBv(1:NewNFREQ,1:NLINES,-1:nAlts+2), &
                   INTEN(1:NLINES,-1:nAlts+2), &
                   NHCNM(1:NLINES,-1:nAlts+2), &
                  gheat1(1:NLINES,-1:nAlts+2), &
             NewTAU(1:NewNFREQ,1:NLINES,-1:nAlts+2), &
             NewVOI(1:NewNFREQ,1:NLINES,-1:nAlts+2), &
                   VULIM(1:NLINES,-1:nAlts+2), &
                   VLLIM(1:NLINES,-1:nAlts+2) )

       PRODHEAT1(1:NLINES,-1:nAlts+2) =  &
                        gheat1(1:NLINES,-1:nAlts+2)* &
                        INTEN(1:NLINES,-1:nAlts+2)*(2.0*PI)
 
       SMHEAT1(-1:nAlts+2) = &
                 SUM(PRODHEAT1(1:NLINES,-1:nAlts+2),1)
 
       HEAT1(-1:nAlts+2)  =  &
                       SMHEAT1(-1:nAlts+2)* &
                          NHCN(-1:nAlts+2)*1.0e+06  ! J/(m^3 s) 

        CALL Newheat4gauss( &
                NewV(1:NewNFREQ,1:NLINES,-1:nAlts+2), &
                    INTEN(1:NLINES,-1:nAlts+2), &
                   ALPHAD(1:NLINES,-1:nAlts+2), &
                    NHCNM(1:NLINES,-1:nAlts+2), &
                   gheat4(1:NLINES,-1:nAlts+2), &
                       TM(1:NLINES,-1:nAlts+2), &
                       ZM(1:NLINES,-1:nAlts+2), &
                      VJM(1:NLINES,-1:nAlts+2), &
              NewTAU(1:NewNFREQ,1:NLINES,-1:nAlts+2), &
              NewVOI(1:NewNFREQ,1:NLINES,-1:nAlts+2), &
                    VULIM(1:NLINES,-1:nAlts+2), &
                    VLLIM(1:NLINES,-1:nAlts+2), &
                    NewBv(1:NewNFREQ,1:NLINES,-1:nAlts+2) )
     !  call cpu_time(timeend)
!     !  print '("New Heat 4 Time=",f6.3,"seconds.")',timeend - timestart
!  
        PRODHEAT4(1:NLINES,-1:nAlts+2) =  &
                         gheat4(1:NLINES,-1:nAlts+2)* &
                         INTEN(1:NLINES,-1:nAlts+2)*(2.0*PI)
 
        SMHEAT4(-1:nAlts+2) = &
                  SUM(PRODHEAT4(1:NLINES,-1:nAlts+2),1)

        HEAT4(-1:nAlts+2)  =  &
                        SMHEAT4(-1:nAlts+2)* &
                           NHCN(-1:nAlts+2)*1.0E+06      ! J/(m^3 s) 

     RadCoolingRate(iLon,iLat,1:nAlts,iBlock) = &
                 COOL(1:nAlts) - &
                HEAT1(1:nAlts) - &
                HEAT4(1:nAlts) 

     HCNCoolRatio(iLon,iLat,1:nAlts,iBlock) = &
            RadCoolingRate(iLon,iLat,1:nAlts,iBlock)/COOL(1:nAlts)

!   real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: HCNCoolRatio

     enddo  ! End Outer Latitude-Loop
   enddo  ! End Outer Longitude-Loop

  call end_timing("calc_radcooling")

      contains

  subroutine Newheat1gauss(Planck,SJ,NHCNM,gheat1,TauIn,PHI,B,A)

      implicit none
!
! INPUT Variables:------------------------------------------
!
      REAL, INTENT(IN), DIMENSION(NLINES,-1:nAlts+2) ::   &
       SJ,      &
       NHCNM,   &
       B,       &
       A  
!
      REAL,INTENT(IN),DIMENSION(NewNFREQ,NLINES,-1:nAlts+2) ::  &
       Planck,  &
       TauIn,  &
       PHI
!
      REAL, INTENT(OUT), DIMENSION(NLINES,-1:nAlts+2) ::   &
       gheat1
!
! LOCAL Variables:------------------------------------------
!
      INTEGER :: i,j,k,m
      INTEGER :: iAlt,iLine,iLon,iFreq
      INTEGER :: iiAlt
!
      REAL, DIMENSION(NewNFREQ,NLINES,-1:nAlts+2) ::   &
       TOTAL

      REAL ::   &
       TOL
!
      REAL, DIMENSION(NLINES,-1:nAlts+2) ::   &
       TOTAL0,  &
       FACTOR20
!
!!!! THESE VARIALBES LOSE ONE ALT INDEX

      REAL, DIMENSION(NLINES,-1:nAlts + 2) ::   &  
       TAUUP,  &
       TAUDWN,  &
       DTau,  &
       TauMax, &
       EXPTAU

      REAL, DIMENSION(NLINES,-1:nAlts+2) ::  &
        S,  &
        OLD 

        TOTAL0 = 0.0
        TOTAL = 0.0
!
! heat1gauss begin:------------------------------------------
      i = 1
      j = 1
      k = 1
      m = 1
      iiAlt = 1
      TOL = 1.0e-5

      FACTOR20       = 0.0
      FACTOR20(:,-1) = 1.0

      do iFreq = 1, NewNFREQ   ! Outer Most Loop over Frequencies

          do iAlt = -1, nAlts + 2
            TauMax(1:NLINES,iAlt) = TauIn(iFreq,1:NLINES,-1)   
          enddo 
          DTAU(:,:) = TauMax(:,:) - TauIn(iFreq,:,:)  
          call expint1(DTau(:,0:nAlts+2),EXPTAU(:,0:nAlts+2),nAlts)
!
          FACTOR20(:,-1) = 1.0
          do iAlt =  0, nAlts + 2
               FACTOR20(:,iAlt) =    &
                EXP(-DTAU(:,iAlt) ) -  &
                     DTAU(:,iAlt)*EXPTAU(:,iAlt) 
          enddo 

          do iAlt = -1, nAlts + 2
                 TOTAL(iFreq,:,iAlt) =  &
                    FACTOR20(:,iAlt)*  &
                    PHI(iFreq,:,iAlt)
          enddo 
      enddo ! Outer iFreq loop

      S = 0.0   ! This is our Integral Sum
      S = S + FreqQuadWeights(1)*&
          (Planck(MiddleFreqIndex,:,:)*TOTAL(MiddleFreqIndex,:,:))
      ! We assume that the line-shape is symmetric in Frequency

      do iFreq = 1, FreqIndexMax
          S = S + FreqQuadWeights(iFreq + 1)*&
             (Planck(MiddleFreqIndex+iFreq,:,:)*TOTAL(MiddleFreqIndex+iFreq,:,:) + &
              Planck(MiddleFreqIndex-iFreq,:,:)*TOTAL(MiddleFreqIndex-iFreq,:,:))
      enddo 

      gheat1 = (0.5)*(B-A)*S

  end subroutine Newheat1gauss


! --------------------------------------------------------------------------
! Sumfgausstau declaration
!---------------------------------------------------------------------------

      subroutine Sumfgausstau(SJ,NHCNM,Z,PHI,ITAU,H_HCN)

      IMPLICIT NONE

!
! INPUTS ----------------------------------------------
!
      REAL,INTENT(IN),DIMENSION(-1:nAlts+2) ::  &
       Z
!
      REAL,INTENT(IN),DIMENSION(1:NLINES,-1:nAlts+2) ::  &
       SJ,  &
       NHCNM
!
      REAL,INTENT(IN),DIMENSION(1:NewNFREQ,1:NLINES,-1:nAlts+2) ::  &
       PHI

      REAL,INTENT(OUT),DIMENSION(1:NewNFREQ,1:NLINES,-1:nAlts+2)::  &
       ITAU 

      real, intent(in) :: H_HCN

! END INPUTS ----------------------------------------------
!
! Local Vars ----------------------------------------------
      INTEGER ::   &
       j,  &
       k,  &
       m,  &
       n,  &
       i,  &
       loc

      INTEGER ::   &
       iAlt, &
       iLon, &
       iLine

      INTEGER ::   &
       istart,  &
       iend,  &
       irate

      REAL,DIMENSION(NLINES,-1:nAlts+2) ::  &
       TOTAL, LnTotal

      REAL,DIMENSION(NLINES) ::  &
       NewTOTALAbove, NewLnTotalAbove, &
       NewTOTALBelow, NewLnTotalBelow

      REAL,DIMENSION(NLINES) ::  &
       InterpLnTotal
!
      REAL,DIMENSION(NLINES,-1:nAlts+2) ::  &
       ITERP2 

      REAL,DIMENSION(NPTS,NewNFREQ,NLINES,-1:nAlts+2) ::  &
       DEPTH
!
      real,dimension(NPTS) :: x
!
      real :: A, B
      real :: TopDepth(NewNFREQ,NLINES)   !!! TAU At the top
      real :: Old(1:NLINES), TauSum(1:NLINES)
!
!----------END Local Vars-----------------------
      m = 1
      DO 
      IF (2**m .EQ. NPTS) EXIT
         m = m + 1
      END DO
!      write(*,*) 'm = ', m
!!!! We integrate from the top downward
!!!! We can then add each increment to the next
!!!! to get the total
      ITAU = 0.0
      do iFreq = 1, NewNFREQ
        do iLine = 1, NLINES
           TopDepth(iFreq,iLine) = &   !!! TAU At the top
                PHI(iFreq,iLine,nAlts+2)*&   ! 1/Hz
                       SJ(iLine,nAlts+2)*&   ! Hz*m^2/molecules
                    NHCNM(iLine,nAlts+2)*(1.0e+06)*& ! Convert NHCN from (cm^-3) -> m^-3
                    H_HCN  ! in m
!!! Store TopDepth in our Optical Depth Array
            ITAU(iFreq,iLine,nAlts+2) = TopDepth(iFreq,iLine)  !! Top Optical Depth
!            write(*,*) 'iFreq, iLine, TopDepth = ', iFreq, iLine, TopDepth(iFreq,iLine)
        enddo  ! Rotational Lines Loop
      enddo  ! Frequency Gaussian Points Loop 


!!!! Begin Altitude Loop -------- Integrate Downward
     do iAlt = nAlts+1, -1, -1
        A = Z(iAlt  )*1000.0   !!! Convert to m
        B = Z(iAlt+1)*1000.0   !!! Convert to m
! Now Our Endpoints are B, A which are NLVLS in Size
!        write(*,*) 'iAlt, Bottom(A), and Top (B) = ',&
!                    iAlt, A, B
! These Quadruature Points do not depend upon 
! Frequency or Line number
! They only depend upon the physical space we sample
        do k = 1 , NPTS/2
           x(2*k - 1)= 0.5*( Qf1(k,m)*(B - A) + B + A)
           x(2*k    )= 0.5*(-Qf1(k,m)*(B - A) + B + A)
        enddo 

        ! Do the Middle Frequency (the core of the lineshape)
        NewTotalAbove(1:NLINES) = PHI(MiddleFreqIndex,1:NLINES,iAlt+1)*&
                                                   SJ(1:NLINES,iAlt+1)*&
                                                NHCNM(1:NLINES,iAlt+1)*1.0e+06
        NewLnTotalAbove(1:NLINES) = alog(NewTotalAbove)

        NewTotalBelow(1:NLINES) = PHI(MiddleFreqIndex,1:NLINES,iAlt)*&
                                                   SJ(1:NLINES,iAlt)*&
                                                NHCNM(1:NLINES,iAlt)*1.0e+06
        NewLnTotalBelow(1:NLINES) = alog(NewTotalBelow)


!          write(*,*) 'Testing the Interpolation for Gaussian Points'

! Factors of 1000.0 convert km -> m
        do n = 1, NPTS
           InterpLnTotal(1:NLINES) = NewLnTotalBelow(1:NLINES) + &
             (x(n) - Z(iAlt)*1000.0)* &
             ((NewLnTotalAbove(:) - NewLnTotalBelow(:))/&
             (Z(iAlt+1)*1000.0 - Z(iAlt)*1000.0))  

           DEPTH(n,MiddleFreqIndex,1:NLINES,iAlt) = &
                 exp(InterpLnTotal(1:NLINES))
        enddo 

!!! Optical Depth of the Line Center
        TauSum(1:NLINES) = 0.0
        do k = 1, NPTS/2
          TauSum(1:NLINES) = &
             TauSum(1:NLINES) + &
             Qd(k,m)*(DEPTH(2*k - 1,MiddleFreqIndex,1:NLINES,iAlt) +  &
                      DEPTH(2*k    ,MiddleFreqIndex,1:NLINES,iAlt))
        enddo ! Loop over Altitude Points
            TauSum(1:NLINES) = 0.5*TauSum(1:NLINES)*(B - A)
      
        ITAU(MiddleFreqIndex,1:NLINES,iAlt  ) = &
        ITAU(MiddleFreqIndex,1:NLINES,iAlt+1) + &
                      TauSum(1:NLINES)

!!! Iterate over the different Frequencies
!!! We assume symmetry about the central Index
!!! Due to the symmetry in the lineshapes.
!!! The Central Index is MiddleFreqIndex

        do iFreq = 1, FreqIndexMax
           NewTotalAbove(1:NLINES) = PHI(MiddleFreqIndex+iFreq,1:NLINES,iAlt+1)*&
                                     SJ(1:NLINES,iAlt+1)*&
                                  NHCNM(1:NLINES,iAlt+1)*1.0e+06
           NewLnTotalAbove(1:NLINES) = alog(NewTotalAbove)

           NewTotalBelow(1:NLINES) = PHI(MiddleFreqIndex+iFreq,1:NLINES,iAlt)*&
                                     SJ(1:NLINES,iAlt)*&
                                  NHCNM(1:NLINES,iAlt)*1.0e+06
           NewLnTotalBelow(1:NLINES) = alog(NewTotalBelow)

           do n = 1, NPTS
              InterpLnTotal(1:NLINES) = NewLnTotalBelow(1:NLINES) + &
                (x(n) - Z(iAlt)*1000.0)* &
                ((NewLnTotalAbove(:) - NewLnTotalBelow(:))/&
                (Z(iAlt+1)*1000.0 - Z(iAlt)*1000.0))  

!!! The Function is Symmetric about the Middle Frequency
              DEPTH(n,MiddleFreqIndex + iFreq,1:NLINES,iAlt) = &
                exp(InterpLnTotal(1:NLINES))

!! Now we have our Depth Gaussian point
              DEPTH(n,MiddleFreqIndex - iFreq,1:NLINES,iAlt) = &
                exp(InterpLnTotal(1:NLINES))

           enddo  ! Loop over Altitude Gaussian Points

!!! Now Calculate the Optical Depth of the new layer (iAlt)
!!! And add the layer above it (iAlt + 1)
           TauSum(1:NLINES) = 0.0
           do k = 1, NPTS/2
             TauSum(1:NLINES) = &
                TauSum(1:NLINES) + &
                Qd(k,m)*(DEPTH(2*k - 1,MiddleFreqIndex + iFreq,1:NLINES,iAlt) +  &
                         DEPTH(2*k    ,MiddleFreqIndex + iFreq,1:NLINES,iAlt))
           enddo ! Loop over Altitude Points
               TauSum(1:NLINES) = 0.5*TauSum(1:NLINES)*(B - A)
         
           ITAU(MiddleFreqIndex + iFreq,1:NLINES,iAlt  ) = &
           ITAU(MiddleFreqIndex + iFreq,1:NLINES,iAlt+1) + &
                  TauSum(1:NLINES)

           ITAU(MiddleFreqIndex - iFreq,1:NLINES,iAlt  ) = &
           ITAU(MiddleFreqIndex - iFreq,1:NLINES,iAlt+1) + &
                  TauSum(1:NLINES)
        enddo !k = 1, iFreq, FreqIndexMax
     enddo !iAlt = nAlts+1, -1, -1
    end subroutine Sumfgausstau

      subroutine Newheat4gauss(V,SJ,ALPHAD,NHCNM,gheat4,TM,  &
                           ZM,VJM,TAU,PHI,VULIM,VLLIM, &
                           Planck)

      IMPLICIT NONE
!
! INPUT Variables:------------------------------------------
!
      REAL, INTENT(IN), DIMENSION(NLINES,-1:nAlts+2) ::  &
       ALPHAD,  &
       SJ,  &
       NHCNM,  &
       TM,  &
       ZM,  &
       VJM,  &
       VULIM,  &
       VLLIM  
!
      REAL,INTENT(IN),DIMENSION(NewNFREQ,NLINES,-1:nAlts+2) ::  &
       V,  &
       TAU,  &
       Planck,  &
       PHI
!
      REAL, INTENT(OUT), DIMENSION(NLINES,-1:nAlts+2) ::  &
       gheat4
!
! LOCAL Variables:------------------------------------------

      INTEGER :: i,j,k,m,p,t,&
                 loc,ilon,iline

      INTEGER ::  & 
       start,  &
       end,  &
       rate 

      REAL,DIMENSION(NLINES,-1:nAlts+2) ::  &
       ATAU,  &
       BTAU,  &
       A1TAU,  &
       B1TAU

      REAL,DIMENSION(NLINES,-1:nAlts+2) ::   &
       ALPHA,   &
       B,   &
       FACTOR,   &
       BH,   &
       BL,   &
       SH,   &
       SL,   &
       OLDH,   &
       OLDL,   &
       S,   &
       OLD
!
      REAL,DIMENSION(NPTS,NLINES,-1:nAlts+2) ::  &
       XH,  & ! Higher Integral (above the value of Z)
       XL,  &  ! Lower Integral (below the value of Z)
       TOTH,  & ! Higher Integral (above the value of Z)
       TOTL,  &  ! Lower Integral (below the value of Z)
       TESTH,  &
       TESTL,  &
       TAUV,  &
       EXPH,  &
       EXPL
!
      REAL,DIMENSION(NewNFREQ,NLINES,-1:nAlts+2) ::  &
       QTOT 

      REAL,DIMENSION(NewNFREQ,NLINES,-1:nAlts+2) ::  &
       ContFuncBelow, &  
       ContFuncAbove, &
       TauMax,        &
       LnTau,         &
       LnBv,         &
       IntTauBelow,   &
       IntTauAbove  

      REAL,DIMENSION(NPTS) ::  &
       LocalXH

      REAL,DIMENSION(NLINES) ::  &
       InterpLnTau, &
       InterpLnBv,  &
       InterpTau, &
       DInterpTau, &
       E1Tau, &
       InterpBv

      REAL,DIMENSION(-1:nAlts+2) ::  &
       Zloc

      REAL,DIMENSION(NPTS,NLINES) ::  &
       Integrand

      integer :: iAlt, jAlt, kAlt

      iline = 1
      i = 1
      Zloc = ZM(35,-1:nAlts+2)
!      write(*,*) 'Zloc = ', Zloc

!! Pick out the Proper Quadrature Points
      m = 1
      DO 
        IF (2**m .EQ. NPTS) EXIT
        m = m + 1
      END DO

! ---------------
! FOR LOOP BEGIN-
! ---------------
!!! Need two integrations for this term
!!! One over altitude--(inner integration)
!!! Outer one over frequency

   LnTau = alog(TAU)
    LnBv = alog(Planck)
   ContFuncAbove(:,:,:) = 0.0
   ContFuncBelow(:,:,:) = 0.0

   ! First, Grab the Core Frequency

   do k = 0, FreqIndexMax
      iFreq = MiddleFreqIndex + k

      ! Outer Altitude Loop (the Z in our Integral)
      ! We don't need to iterate over -1 or nAlts + 2
      ! Just go from 0:nAlts+1
      do iAlt =  1, nAlts 
!         write(*,*) 'Outer Alt Loop = ', iAlt, Zloc(iAlt)
         Integrand(:,:) = 0.0
         ContFuncBelow(MiddleFreqIndex + k,:,iAlt) = 0.0
         ContFuncBelow(MiddleFreqIndex - k,:,iAlt) = 0.0
         ContFuncAbove(MiddleFreqIndex + k,:,iAlt) = 0.0
         ContFuncAbove(MiddleFreqIndex - k,:,iAlt) = 0.0

         do jAlt = -1, iAlt-1   ! This is our Z''
            !if (iAlt .eq. -1) exit 
            ! We don't execute this at all for iAlt = -1 
            ! Gives the Gaussain Points in terms of TAU Values

            !write(*,*) 'Contribution Below You, jAlt = ', jAlt
            !write(*,*) 'Zloc(jAlt), Zloc(iAlt-1)= ', &
            !            Zloc(jAlt), Zloc(jAlt+1)
            do p = 1 , NPTS/2
               ! Get integration points in the slab of atmosphere below us.
               ! Should be (NPTS,1:NLINES) in Size 
               ! Exploit the Symmetry of the X(2*p) and X(2*p - 1)
               LocalXH(2*p - 1) = &
                 0.5*( Qf1(p,m)*(Zloc(jAlt+1) - Zloc(jAlt))&
                   + Zloc(jAlt) + Zloc(jAlt+1) )

               LocalXH(2*p    ) = &
                 0.5*(-Qf1(p,m)*(Zloc(jAlt+1) - Zloc(jAlt))&
                   + Zloc(jAlt) + Zloc(jAlt+1) )
               !write(*,*) 'GaussPoint(2p), GaussPoint(2p-1)  = ', &
               !           p, LocalXH(2*p-1), LocalXH(2*p)

            enddo
          ! At each Gausspoint, we need to do a linear interpolation:
!              write(*,*) 'Tau(jAlt) and Tau(jAlt+1) = ', &
!                  TAU(iFreq,35,jAlt), TAU(iFreq,35,jAlt+1)

!             write(*,*) 'Planck(jAlt) and Planck(jAlt+1) = ', &
!                 Planck(iFreq,35,jAlt), Planck(iFreq,35,jAlt+1)

            do p = 1 , NPTS
               InterpLnTau = LnTau(iFreq,1:NLINES,jAlt) + &
                        (LocalXH(p) - Zloc(jAlt))*&
                        (LnTau(iFreq,1:NLINES,jAlt+1) - LnTau(iFreq,1:NLINES,jAlt))/&
                        (Zloc(jAlt+1) - Zloc(jAlt))
               InterpTau = exp(InterpLnTau)   ! 1:NLINES Long

               DInterpTau = InterpTau(1:NLINES) - TAU(iFreq,1:NLINES,iAlt) 


               InterpLnBv = LnBv(iFreq,1:NLINES,jAlt) + &
                        (LocalXH(p) - Zloc(jAlt))*&
                        (LnBv(iFreq,1:NLINES,jAlt+1) - LnBv(iFreq,1:NLINES,jAlt))/&
                        (Zloc(jAlt+1) - Zloc(jAlt))

               InterpBv = exp(InterpLnBv)   ! 1:NLINES Long

!              write(*,*) 'p, InterpBv = ', p, InterpBv(35)


!              write(*,*) 'p, DInterpTau = ', p, DInterpTau(35)
               CALL expint4(DInterpTau(:),E1Tau(:))
!              write(*,*) 'p, DTau, ExpInt(DTau) = ',p,DInterpTau(35),E1Tau(35) 

!               Integrand(p,1:NLINES) = &
!                    E1Tau(:)*InterpBv(:)*(&
!                    TAU(iFreq,1:NLINES,jAlt) - TAU(iFreq,1:NLINES,jAlt+1))

               Integrand(p,1:NLINES) = &
                    E1Tau(:)*InterpBv(:)

               ! The Integrand is Bv*E1(Tau(z) - Tau(z'))*dTau
               ! The dTAu is the size of the "slab" of atmosphere between
               ! jAlt + 1 and jAlt
            enddo 

            !! Add the Altitudes together to get a total contribution from the slab
            ! We need to keep adding each slab to this ConFuncBelow Variable
            do p = 1 , NPTS/2
               ContFuncBelow(iFreq,:,iAlt) =  &
                  ContFuncBelow(iFreq,:,iAlt) + &
                  0.5*Qd(p,m)*(Integrand(2*p  ,1:NLINES) + &
                           Integrand(2*p-1,1:NLINES))*&
                    (TAU(iFreq,1:NLINES,jAlt) - TAU(iFreq,1:NLINES,jAlt+1))
            enddo ! p loop

               ! Assume Symmetry around the Central Frequency
               ContFuncBelow(MiddleFreqIndex - k,:,iAlt) =  &
               ContFuncBelow(iFreq,:,iAlt) 
            ! This Loop integrates over the Interpolation points, 
            ! Which completely accounts for the atmospheric slab
            ! Between jAlt and jAlt+1
        enddo ! Loop over jAlt (Z') Bell et al. [2008] 


      ! Same as jAlt Loop above, but now for the 
      ! Atmosphere above you
        Integrand(:,:) = 0.0
        ContFuncAbove(iFreq,:,iAlt) = 0.0
        do jAlt = iAlt+1, nAlts+2 
!           write(*,*) 'Contribution Above You, jAlt = ', jAlt
!           write(*,*) 'Zloc(jAlt-1), Zloc(iAlt)= ', &
!                       Zloc(jAlt-1), Zloc(jAlt)
           do p = 1 , NPTS/2
              ! Get integration points in the slab of atmosphere above us.
              ! Should be (NPTS,1:NLINES) in Size 
              ! Exploit the Symmetry of the X(2*p) and X(2*p - 1)
              LocalXH(2*p - 1) = &
                0.5*( Qf1(p,m)*(Zloc(jAlt) - Zloc(jAlt-1))&
                  + Zloc(jAlt-1) + Zloc(jAlt) )

              LocalXH(2*p    ) = &
                0.5*(-Qf1(p,m)*(Zloc(jAlt) - Zloc(jAlt-1))&
                  + Zloc(jAlt-1) + Zloc(jAlt) )

!              write(*,*) 'GaussPoint(2p), GaussPoint(2p-1)  = ', &
!                         p, LocalXH(2*p-1), LocalXH(2*p)

           enddo

         ! At each Gausspoint, we need to do a linear interpolation:

!            write(*,*) 'Tau(jAlt-1) and Tau(jAlt) = ', &
!                TAU(iFreq,35,jAlt-1), TAU(iFreq,35,jAlt)

!           write(*,*) 'Planck(jAlt-1) and Planck(jAlt) = ', &
!               Planck(iFreq,35,jAlt-1), Planck(iFreq,35,jAlt)

           do p = 1 , NPTS
              InterpLnTau = LnTau(iFreq,1:NLINES,jAlt-1) + &
                       (LocalXH(p) - Zloc(jAlt-1))*&
                       (LnTau(iFreq,1:NLINES,jAlt) - LnTau(iFreq,1:NLINES,jAlt-1))/&
                       (Zloc(jAlt) - Zloc(jAlt-1))
              InterpTau = exp(InterpLnTau)   ! 1:NLINES Long


 !             write(*,*) 'p, InterpTau = ', p, InterpTau(35)
              DInterpTau = -1.0*(InterpTau(1:NLINES) - TAU(iFreq,1:NLINES,iAlt)) 

              InterpLnBv = LnBv(iFreq,1:NLINES,jAlt-1) + &
                       (LocalXH(p) - Zloc(jAlt-1))*&
                       (LnBv(iFreq,1:NLINES,jAlt) - LnBv(iFreq,1:NLINES,jAlt-1))/&
                       (Zloc(jAlt) - Zloc(jAlt-1))
              InterpBv = exp(InterpLnBv)   ! 1:NLINES Long

!             write(*,*) 'p, DInterpTau = ', p, DInterpTau(35)
             CALL expint4(DInterpTau(:),E1Tau(:))
!             write(*,*) 'p, DTau, ExpInt(DTau) = ',p,DInterpTau(35),E1Tau(35) 

!             Integrand(p,1:NLINES) = &
!                  E1Tau(:)*InterpBv(:)*(&
!                  TAU(iFreq,1:NLINES,jAlt-1) - TAU(iFreq,1:NLINES,jAlt))

             Integrand(p,1:NLINES) = &
                  E1Tau(:)*InterpBv(:)

              ! The Integrand is Bv*E1(Tau(z) - Tau(z'))*dTau
              ! The dTAu is the size of the "slab" of atmosphere between
              ! jAlt + 1 and jAlt
           enddo 

           !! Add the Altitudes together to get a total contribution from the slab
           ! We need to keep adding each slab to this ConFuncBelow Variable
           do p = 1 , NPTS/2
              ContFuncAbove(iFreq,:,iAlt) =  &
                 ContFuncAbove(iFreq,:,iAlt) + &
                 0.5*Qd(p,m)*(Integrand(2*p  ,1:NLINES) + &
                          Integrand(2*p-1,1:NLINES))*&
                ( TAU(iFreq,1:NLINES,jAlt-1) - TAU(iFreq,1:NLINES,jAlt))
     
           enddo ! p loop
               ContFuncAbove(MiddleFreqIndex - k,:,iAlt) =  &
               ContFuncAbove(iFreq,:,iAlt) 
           ! This Loop integrates over the Interpolation points, 
           ! Which completely accounts for the atmospheric slab
           ! Between jAlt and jAlt+1
       enddo ! Loop over jAlt (Z') Bell et al. [2008] 

      enddo !iAlt  Outer Integration Loop
   enddo !k = 0,FreqIndexMax

   ! NewNFreqs, NLines, -1:nAlts+2
   QTOT(:,:,1:nAlts) = &
     ContFuncBelow(:,:,1:nAlts) + ContFuncAbove(:,:,1:nAlts)

   QTOT(:,:, 0) = QTOT(:,:,1) 
   QTOT(:,:,-1) = QTOT(:,:,1) 
   QTOT(:,:,nAlts+1) = QTOT(:,:,nAlts) 
   QTOT(:,:,nAlts+2) = QTOT(:,:,nAlts) 

   S = 0.0
   S = S + FreqQuadWeights(1)*(QTOT(MiddleFreqIndex,:,:)*PHI(MiddleFreqIndex,:,:) )

   do iFreq = 1, FreqIndexMax
      S = S + FreqQuadWeights(iFreq + 1)*&
               (QTOT(MiddleFreqIndex + iFreq,:,:)*PHI(MiddleFreqIndex + iFreq,:,:) + &
                QTOT(MiddleFreqIndex - iFreq,:,:)*PHI(MiddleFreqIndex - iFreq,:,:) )
   enddo 
   gheat4 = 0.5*(VULIM-VLLIM)*S

   end subroutine Newheat4gauss


      subroutine expint1(x,ans,NLVLS)
      implicit none

      INTEGER, INTENT(IN) :: NLVLS 
      REAL, INTENT(IN ), DIMENSION(NLINES, 0:NLVLS+2) ::  x
      REAL, INTENT(OUT), DIMENSION(NLINES, 0:NLVLS+2) :: ans
      
      INTEGER, PARAMETER :: MAXIT = 100 
!      INTEGER, PARAMETER :: MAXIT = 10
!      INTEGER, PARAMETER :: MAXIT = 5 
!      INTEGER, PARAMETER :: MAXIT = 2 
!      INTEGER, PARAMETER :: MAXIT = 1 
      INTEGER :: i,j,k 

      REAL, PARAMETER :: EULER = .577215664901532860D0 
      REAL, PARAMETER :: EPS = epsilon(x(1,1)) 
      REAL, PARAMETER :: BIG = huge(x(1,1))*EPS 

      REAL, DIMENSION(NLINES, 0:NLVLS+2) ::  ID

      REAL, DIMENSION(NLINES, 0:NLVLS+2) ::   &
       a,  &
       b,  &
       c,  &
       d,  &
       h,  &
       del,  &
       fact,  &
       EXPH,  &
       EXPL

       ID = 1.0
! ---------------------------------------------------------
! BEGIN CALCULATION:---------------------------------------
! FORM FOR X >= 1.0 
! ---------------------------------------------------------
       b = x + ID
       c = BIG*ID 
       d = ID/b 
       h = d 
       DO i = 1,MAXIT
          a = -1.0D0*i*i*ID
          b = b + 2.0D0*ID
          d = ID/(a*d + b) 
          c = b + a/c
          del = c*d
          h = h*del
          IF ( MAXVAL(ABS(del) - 1.0) <= EPS*1.0D13) EXIT
       ENDDO

       EXPH = h*exp(-x)

! ---------------------------------------------------------
! BEGIN CALCULATION:---------------------------------------
! FORM FOR X < 1.0 
! ---------------------------------------------------------

       EXPL = -alog(x) - EULER 

       fact = 1.0D0 
       del = 0.0D0 
       DO i = 1, MAXIT
          fact = -(fact*x)/i 
          del = -fact/i
          EXPL = EXPL + del
          IF ( MAXVAL(ABS(del)) <= MAXVAL(ABS(EXPL))*EPS) EXIT
       ENDDO

       WHERE(x < 1.0D0)
         ans = EXPL
       ELSEWHERE
         ans = EXPH
       END WHERE
      END subroutine expint1


      subroutine expint4(x,ans)
      implicit none

      REAL, INTENT(IN), DIMENSION(NLINES) ::  x
      REAL, INTENT(OUT), DIMENSION(NLINES) :: ans
      

!      INTEGER, PARAMETER :: MAXIT = 1 
!      INTEGER, PARAMETER :: MAXIT = 10
      INTEGER, PARAMETER :: MAXIT = 5 
!     INTEGER, PARAMETER :: MAXIT = 2 
!      INTEGER, PARAMETER :: MAXIT = 1 
      INTEGER :: i,j,k 

      REAL, PARAMETER :: EULER = .577215664901532860D0 
      REAL, PARAMETER :: EPS = epsilon(x(1)) 
      REAL, PARAMETER :: BIG = huge(x(1))*EPS 

      REAL, DIMENSION(NLINES) ::  ID

      REAL, DIMENSION(NLINES) ::   &
       a,  &
       b,  &
       c,  &
       d,  &
       h,  &
       del,  &
       fact,  &
       EXPH,  &
       EXPL

       ID = 1.0
! ---------------------------------------------------------
! BEGIN CALCULATION:---------------------------------------
! FORM FOR X >= 1.0 
! ---------------------------------------------------------
       b = x + ID
       c = BIG*ID 
       d = ID/b 
       h = d 
       DO i = 1,MAXIT
          a = -1.0D0*i*i*ID
          b = b + 2.0D0*ID
          d = ID/(a*d + b) 
          c = b + a/c
          del = c*d
          h = h*del
          IF ( MAXVAL(ABS(del) - 1.0) <= EPS*1.0e13) EXIT
       ENDDO

       EXPH = h*exp(-x)

! ---------------------------------------------------------
! BEGIN CALCULATION:---------------------------------------
! FORM FOR X < 1.0 
! ---------------------------------------------------------

       EXPL = -alog(x) - EULER 

       fact = 1.0
        del = 0.0
       DO i = 1, MAXIT
          fact = -(fact*x)/i 
          del = -fact/i
          EXPL = EXPL + del
          IF ( MAXVAL(ABS(del)) <= MAXVAL(ABS(EXPL))*EPS) EXIT
       ENDDO

       WHERE(x < 1.0)
         ans = EXPL
       ELSEWHERE
         ans = EXPH
       END WHERE

      END subroutine expint4

     end subroutine calc_radcooling 




   subroutine calc_saturn_tides(iBlock)
     use ModPlanet
     use ModGITM
     use ModConstants
     use ModTime, only : tSimulation

     integer, intent(in) :: iBlock

     real :: a_titan, e_titan, Mass_Saturn 
     real :: omega_titan
     real :: GravConstant
     real :: GravScalar
     integer :: iAlt, iLat, iLon, iDir

     a_titan = 1.22183e+09   ! Titan orbital semi-major axis (in m) 
     e_titan = 0.0292        ! Titan orbital eccentricity (e)
     omega_titan = OMEGABodyInput ! Titan orbital angular frequency (rads/s)
     Mass_Saturn = 5.685e+26 ! Mass of Saturn
     GravConstant = 6.67259e-11  ! Graviational constant in SI units (N m^2/kg^2)
     !Note:  Titan is phase-locked so that it's rotation rate and orbit rate are the same

     !This Grav Scalar occurs so frequently, we just calculate it once and use it below
     GravScalar = GravConstant*Mass_Saturn/(a_titan**3.0)


     do iAlt = -1, nAlts+2
       do iLat = -1, nLats+2
         do iLon = -1, nLons+2
            ! ----------------
            ! Time-Independent Component of the Tides
            ! This would be the only component if Titan were in a perfectly
            ! Circular Orbit
            ! ----------------
            TideAIndp(iLon,iLat,iAlt,iUp_,iBlock) = GravScalar*&
                  RadialDistance_GB(iLon,iLat,iAlt,iBlock)*&
                 (3.0*cos(Latitude(iLat,iBlock))*cos(Latitude(iLat,iBlock))*&
                      cos(Longitude(iLon,iBlock))*cos(Longitude(iLon,iBlock)) - &
                  1.0 )

            TideAIndp(iLon,iLat,iAlt,iNorth_,iBlock) = -1.0*GravScalar*&
                  RadialDistance_GB(iLon,iLat,iAlt,iBlock)*&
                 (3.0*cos(Latitude(iLat,iBlock))*sin(Latitude(iLat,iBlock))*&
                     cos(Longitude(iLon,iBlock))*cos(Longitude(iLon,iBlock)) )

            TideAIndp(iLon,iLat,iAlt,iEast_,iBlock) = -1.0*GravScalar*&
                  RadialDistance_GB(iLon,iLat,iAlt,iBlock)*&
                 (3.0*cos(Latitude(iLat,iBlock))*&
                     cos(Longitude(iLon,iBlock))*sin(Longitude(iLon,iBlock)) )

            ! Radial Tide Component:  Note proportional to Titan's eccentricity, e
            ! This is due to Titan moving closer and away from Saturn in its orbit
            ! We assume that we start at Periapsis at t = 0.0
            ! By JPL Ephemeris, this occurs on Jan, 04, 2007, 1500 hours
            ! The Local Time at Longitude 0 is roughly midnight 
            TideADep(iLon,iLat,iAlt,iUp_,iBlock) = GravScalar*&
                  RadialDistance_GB(iLon,iLat,iAlt,iBlock)*&
                 (3.0*cos(Latitude(iLat,iBlock))*cos(Latitude(iLat,iBlock))*&
                      cos(Longitude(iLon,iBlock))*cos(Longitude(iLon,iBlock)) - &
                  1.0 )*(3.0*e_titan*cos(OmegaBodyInput*tSimulation))
              
            TideADep(iLon,iLat,iAlt,iNorth_,iBlock) = -1.0*GravScalar*&
                  RadialDistance_GB(iLon,iLat,iAlt,iBlock)*&
                 (3.0*cos(Latitude(iLat,iBlock))*sin(Latitude(iLat,iBlock))*&
                     cos(Longitude(iLon,iBlock))*cos(Longitude(iLon,iBlock)) )*&
                 (3.0*e_titan*cos(OmegaBodyInput*tSimulation))
!
            TideADep(iLon,iLat,iAlt,iEast_,iBlock) = -1.0*GravScalar*&
                  RadialDistance_GB(iLon,iLat,iAlt,iBlock)*&
                 (3.0*cos(Latitude(iLat,iBlock))*&
                     cos(Longitude(iLon,iBlock))*sin(Longitude(iLon,iBlock)) )*&
                 (3.0*e_titan*cos(OmegaBodyInput*tSimulation))
 
!             ! Add the Librational Tide to the Time-dependent part.
!             ! This is due to the motion of Saturn in longitude about the Longitude 
!             ! of zero.
            TideADep(iLon,iLat,iAlt,iUp_,iBlock) = &
            TideADep(iLon,iLat,iAlt,iUp_,iBlock) + 1.0*GravScalar*&
                  RadialDistance_GB(iLon,iLat,iAlt,iBlock)*&
                 (6.0*e_titan*cos(Latitude(iLat,iBlock))*cos(Latitude(iLat,iBlock))*&
                      sin(2.0*Longitude(iLon,iBlock)))*sin(OmegaBodyInput*tSimulation)
!                 
            TideADep(iLon,iLat,iAlt,iNorth_,iBlock) = &
            TideADep(iLon,iLat,iAlt,iNorth_,iBlock) - 1.0*GravScalar*&
                  RadialDistance_GB(iLon,iLat,iAlt,iBlock)*&
                 (6.0*e_titan*cos(Latitude(iLat,iBlock))*sin(Latitude(iLat,iBlock))*&
                     sin(2.0*Longitude(iLon,iBlock)) )*sin(OmegaBodyInput*tSimulation)
!
            TideADep(iLon,iLat,iAlt,iEast_,iBlock) = &
            TideADep(iLon,iLat,iAlt,iEast_,iBlock) + 1.0*GravScalar*&
                  RadialDistance_GB(iLon,iLat,iAlt,iBlock)*&
                 (6.0*e_titan*cos(Latitude(iLat,iBlock))*cos(2.0*Longitude(iLon,iBlock)))*&
                  sin(OmegaBodyInput*tSimulation)

         enddo !iLon = -1, nLons+2
       enddo !iLat = -1, nLats+2
     enddo !iAlt = -1, nAlts+2

   ! Remove the tides for now
   TideADep = 0.0
   TideAInDep = 0.0

   end subroutine calc_saturn_tides


   subroutine calc_cooltospace(iBlock)
!  ---------------------------------------------------------------------
!  Purpose: 
!  =======
!  This is the Cooling Routine to be used for the Titan TGCM
!  This version incorporates all heating and cooling terms into
!  the full line by line radiative transfer calculations.
!     Date            Programmer                   Version
!  ==========     ==================            ==============
!  09-08-04           Jared Bell                 Original(1.0) 
!  12-02-04           Jared Bell                 Update (1.1)   
!  12-07-04           Jared Bell                 Update (1.2)   
!  01-27-05           Jared Bell                 Consistency Check   
!  01-27-05           Jared Bell                 Added Heating 1 
!  01-27-05           Jared Bell                 Added Heating 2 
!  01-28-05           Jared Bell                 Added Dynamic Memory 
!  09-07-05           Jared Bell                 Modify to use TAU 
!  09-12-05           Jared Bell                 Utilized Linear Interp. 
!  09-16-05           Jared Bell                 Expint more accurate 
!                                                and more efficient 
!  03-02-06           Jared Bell                 GITM                 
!  0l-07-07           Jared Bell                 Computational Efficiency
!                                                Update 
!  05-28-11           Jared Bell                 Removed all 5-D Arrays 
!                                                (nLon, nLat, nAlts, nLines, nFreqas)
!                                                Implemented Loops Instead
!  -----------------------------------------------------------------

      use ModPlanet
      use ModGITM
      use ModConstants, only:  Boltzmanns_Constant, Speed_Light, &
                               Planck_Constant
      use ModInputs, only:  iDebugLevel
      use ModSources, only: RadCoolingRate
      use ModTime, only: tSimulation, CurrentTime
      use ModIndicesInterfaces

      implicit none

! INPUTS:----------------------------------------------------------

      integer,intent(in):: iBlock 

! LOCAL VARS:------------------------------------------------------

      integer :: i,j,l,m,k,n,ii,tj 
      integer :: iLat, iLon, iAlt, iLine, iFreq

      integer, parameter :: NLVLS = nAlts+4

      integer, parameter :: NLINES = rotlines 

!      integer, parameter :: NFREQ = rotfreqs 
      integer, parameter :: NewNFREQ = newrotfreqs 
      integer, parameter :: NPTS = rotpts 

!      integer, parameter :: blon = 1 
!      integer, parameter :: elon = nLons 
!      integer, parameter :: blat = 1 
!      integer, parameter :: elat = nLats 
!
!      real,dimension(-1:nAlts+2) :: &
!       TCOOLa

      real ::   &
       FTH

      real,dimension(-1:nAlts+2) ::    &
       T,       & 
       NHCN,    &
       COOL

      real,dimension(1:NLINES,-1:nAlts+2) ::  &
       gcool,   &
       INTEN,   &
       ALPHAD,  &
       ALPHAV,  &
       ALPHAL,  &    ! VOIGT
       DELTAV,  &
       VULIM,   &
       VLLIM,   &
       VJM,     &
       TM,      &
       ZM,      &
       VTH,     &
       NHCNM,   &
       GammaNew,  &   ! VOIGT
       YRATIO         ! VOIGT

!!!! FWHM of the Doppler and Lorentz Profiles
      real,dimension(1:NLINES,-1:nAlts+2) ::  &
        FWHML, FWHMD, FWHMV

!!! Coefficients for Combining Gaussian and Lorentz shapes
!!! DRatio used by Liu et al. [2001] in JQSRT
      real :: DRatio, CoefL, CoefG

      real,dimension(1:NLINES,-1:nAlts+2) ::  &
       NEW, &
       OLD, &
       Y1,  &
       Y2,  &
       Y

      real,dimension(1:NLINES,-1:nAlts+2) ::  &
       PRODCOOL

      real,dimension(-1:nAlts+2) ::    &
       SMCOOL
!
!

!      real,dimension(-1:nAlts+2) ::  &
!       NHCNT 
!   
! --------------BEGIN TAU VARS--------------------------------------- 
! The Following Variables are used solely within the Tau Integration
! 
      real :: DFACT
      real,dimension(1:NewNFREQ,1:NLINES,-1:nAlts+2) ::  &
       NewDOP,  &  ! the DOPPLER Lineshape
       NewLorentz,  &  ! the DOPPLER Lineshape
       NewX,  &
       NewV, &
       NewYV,  &   ! VOIGT
       NewVOI,  &  ! the VOIGT Lineshape
       NewDopFactor

      real,dimension(1:NewNFREQ,1:NLINES,-1:nAlts+2) ::  &
       NewBv,  &
       NewBalpha

      real :: Bfactor

      integer :: MiddleFreqIndex, FreqIndexMax

      MiddleFreqIndex = 5
      FreqIndexMax = 4

  call start_timing("calc_radcooling")

     do iLon = 1, nLons
       do iLat = 1, nLats

       T(-1:nAlts+2) = &
        Temperature(iLon,iLat,-1:nAlts+2,iBlock)*TempUnit(iLon,iLat,-1:nAlts+2)

       NHCN(-1:nAlts+2) = &
        NDensityS(iLon,iLat,-1:nAlts+2,iHCN_,iBlock)*1.0e-06

! ------------------------------------------------------------------------
! INTERPOLATING INTENSITY DATA TO TEMPERATURE GRID -----------------------
! ------------------------------------------------------------------------
             do iAlt = -1, nAlts+2
              do iLine = 1, NLINES
!\
! Speed_Light*100.0:   Converts cm^-1 --> Hz 
! 1.0e-04          : converts Sj from Hz*cm^2/molecule -> Hz*m^2/molecule
                INTEN(iLine,iAlt) =             &
                  ( Speed_Light*100.0)*(1.0e-4)*          &
                  exp( pmat(iLine,1)                         &
                  + pmat(iLine,2)*(T(iAlt)   )  &
                  + pmat(iLine,3)*(T(iAlt)**2)  &
                  + pmat(iLine,4)*(T(iAlt)**3)  &
                  + pmat(iLine,5)*(T(iAlt)**4)  &
                  + pmat(iLine,6)*(T(iAlt)**5)  &
                  + pmat(iLine,7)*(T(iAlt)**6)  &
                  + pmat(iLine,8)*(T(iAlt)**7)  & 
                  + pmat(iLine,9)*(T(iAlt)**8)  )    


              enddo! iLine = 1, NLINES
             enddo !iAlt = -1, nAlts+2

! ------------------------------------------------------------------------
! Frequency Integration Variables TM,VJM,ZM,NHCNM (1:nlines,1:nAlts)
! ------------------------------------------------------------------------
             do iAlt = -1, nAlts+2 
!              write(*,*) '=====> iAlt =', iAlt
              do iLine = 1, NLINES  

               VJM(iLine,iAlt) = freqhz(iLine)        ! Lines in Hz
                ZM(iLine,iAlt) = Altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0   ! m -> km
                TM(iLine,iAlt) = T(iAlt)       ! Temp in  K
            NHCNM(iLine,iAlt) = NDensityS(iLon,iLat,iAlt,iHCN_,iBlock)*1.0E-06  ! m^-3 -> cm^-3

! These provide the Alpha_L in cm^-1 (we call this GammaNew)
          GammaNew(iLine,iAlt) = &
                    ALPHAL0(iLine)*&
             (1.0 - NDensityS(iLon,iLat,iAlt,iHCN_,iBlock)/NDensity(iLon,iLat,iAlt,iBlock))*&
             (Pressure(iLon,iLat,iAlt,iBlock)/(101325.00))*(298.0/T(iAlt))**ALPHALExp(iLine) + &
                    ALPHAS0(iLine)*&
             (NDensityS(iLon,iLat,iAlt,iHCN_,iBlock)/NDensity(iLon,iLat,iAlt,iBlock))*&
             (Pressure(iLon,iLat,iAlt,iBlock)/(101325.00))*(298.0/T(iAlt))**ALPHALExp(iLine)

!! AlphaL is the Half-width of the Lorentz profile (pressure-broadening) in Hz
          ALPHAL(iLine,iAlt) = &
                 GammaNew(iLine,iAlt)*(Speed_Light*100.0)   !! This is AlphaL in Hz

              enddo  ! End iLine
             enddo ! End iAlt

! ------------------------------------------------------------------------
! DOPPLER HALFWIDTHS AT GIVEN TEMPERATURES -------------------------------
! ------------------------------------------------------------------------
           FTH = SQRT(2.0*(Boltzmanns_Constant/Mass(iHCN_)))
           VTH = FTH*SQRT(TM)
           ALPHAD = (VJM*(VTH/(Speed_Light)))  !!! Our Doppler (Gaussian) half-width
           DELTAV = 5.0*ALPHAD
           VULIM = VJM + DELTAV
           VLLIM = VJM - DELTAV

! =================================================>
! Calculate the FWHM of the Convolved Voigt Profile
! Use the form taken from Olivier and Longbothom

           FWHMD = ALPHAD*(2.0*sqrt(2.0*alog(2.0)))  !! Our FUllWidth at Half-Max of Doppler
           FWHML = 2.0*ALPHAL

           FWHMV = 0.5346*FWHML + &
              sqrt(0.2166*FWHML**2.0 + FWHMD**2.0)

           ALPHAV = FWHMV/2.0  ! Our Voigt Half-width

!----------------
!WHILE LOOP BEGIN
!----------------
!
!      m = 1
!      do 
!      if (2**m .EQ. NFREQ) exit
!         m = m + 1
!      enddo
!----------------
!END WHILE LOOP 
!---------------

         NewV(MiddleFreqIndex,:,:) = &
                  (FreqQuadAbscissas(1)*(VULIM - VLLIM) + VULIM + VLLIM)*.5

       do k = 1 , FreqIndexMax
         NewV(MiddleFreqIndex - k,:,:)= &
              (FreqQuadAbscissas(k+1)*(VULIM - VLLIM) + VULIM + VLLIM)*.5
         NewV(MiddleFreqIndex + k,:,:)=&
              (-1.0*FreqQuadAbscissas(k+1)*(VULIM - VLLIM) + VULIM + VLLIM)*.5
       enddo

         DFACT = SQRT(alog(2.0)/PI)

        do iAlt = -1, nAlts+2
           do iLine = 1, NLINES
              ! Note that DRatio will vary 
              ! from -1.0 ( Pure Doppler)
              ! to 1.0 (Pure Lorentz)
              DRatio = (FWHML(iLine,iAlt) - FWHMD(iLine,iAlt))/&
                       (FWHML(iLine,iAlt) + FWHMD(iLine,iAlt))
 
              ! Coef L and CoefD determine the "relative mixing"
              ! Of our line shapes into the combined Voigt Shape
              CoefL = 0.68188 + 0.61293*DRatio - &
                      0.18384*DRatio**2.0 - &
                      0.11568*DRatio**3.0 

              CoefG = 0.32460 - 0.61825*DRatio + &
                      0.17681*DRatio**2.0 + &
                      0.12109*DRatio**3.0 

              do iFreq = 1, NewNFREQ

              ! Calculate the Doppler contribution
              
              NewX(iFreq,iLine,iAlt) =  &
                  SQRT( alog(2.0) )*  (  ( NewV(iFreq,iLine,iAlt)  - VJM(iLine,iAlt) )/ & 
                   ALPHAV(iLine,iAlt) ) 

              NewDOP(iFreq,iLine,iAlt) =  &
                    (DFACT/ALPHAV(iLine,iAlt))*    & 
                    EXP( -( NewX(iFreq,iLine,iAlt)**2.0) ) 

              NewYV(iFreq,iLine,iAlt) =  &
                    ( NewV(iFreq,iLine,iAlt)  - VJM(iLine,iAlt) )

              NewLorentz(iFreq,iLine,iAlt) =  &
                    (1.0/PI)*ALPHAV(iLine,iAlt)/    & 
                    ( NewYV(iFreq,iLine,iAlt)**2.0 + ALPHAV(iLine,iAlt)**2.0)

              NewVOI(iFreq,iLine,iAlt) = &
                   CoefL*NewLorentz(iFreq,iLine,iAlt) + &
                   CoefG*NewDOP(iFreq,iLine,iAlt) 

 !             write(*,*) '===============> iFreq', iFreq, NewVOI(iFreq,iLine,iAlt) 

              NewBalpha(iFreq,iLine,iAlt)  =             &
                    (  Planck_Constant/Boltzmanns_Constant)*  &
                    (  NewV(iFreq,iLine,iAlt) /        &
                            TM(iLine,iAlt)  )

              Bfactor = 1.0/ &
                        ( EXP(NewBalpha(iFreq,iLine,iAlt) ) - 1.0 ) 

              NewBv(iFreq,iLine,iAlt) =                   &
                     ( ( 2.0*Planck_Constant ) /        & 
                     ( (Speed_Light)**2 )  ) *                   &
                     ( NewV(iFreq,iLine,iAlt)**3 )*Bfactor

         enddo
       enddo
    enddo

! ------------------------------------------------------------------------
! Cooling CALCULATIONS
! Cool-To-Space Portion: 
! ------------------------------------------------------------------------

       OLD = 0.0
       NEW = 0.0
       NEW = NEW + FreqQuadWeights(1)*(NewBv(MiddleFreqIndex, :, :)*NewVOI(MiddleFreqIndex, :, :))

      do iFreq = 1 , FreqIndexMax
        Y1 = NewBv(MiddleFreqIndex - iFreq,:,:)*NewVOI(MiddleFreqIndex - iFreq,:,:)
        Y2 = NewBv(MiddleFreqIndex + iFreq,:,:)*NewVOI(MiddleFreqIndex + iFreq,:,:) 
        NEW = NEW + FreqQuadWeights(iFreq + 1)*(Y1 + Y2)
      enddo

      gcool(:,:) = (VULIM(:,:) - VLLIM(:,:) )*NEW*0.5
! Multiply by the Intensity of the lines and a geometric factor of 4pi
!
      PRODCOOL(1:NLINES,-1:nAlts+2) =  &
                       gcool(1:NLINES,-1:nAlts+2)* &
                       INTEN(1:NLINES,-1:nAlts+2)*(4.0D0*PI)
! Sum over the lines
      SMCOOL(-1:nAlts+2) = &
                SUM(PRODCOOL(1:NLINES,-1:nAlts+2),1)
! Multiply by HCN density 
      COOL(-1:nAlts+2)  =  &
                      SMCOOL(-1:nAlts+2)* &
                        NHCN(-1:nAlts+2)*1.0e+06   ! J/(m^3 s) 

       RadCoolingRate(iLon,iLat,1:nAlts,iBlock) = &
                   COOL(1:nAlts)*HCNCoolRatio(iLon,iLat,1:nAlts,iBlock)

       enddo ! iLat
    enddo ! iLon

   endsubroutine !calc_cooltospace(iBlock)

