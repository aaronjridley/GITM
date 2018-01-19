Module ModChemistry

!  Turn off CO2 sources and losses; infinite reservoir assumption.
!  Fix added for: (1)  O+O2+CO2 and (2) few rate coefficients: 11-12-01
!  Corrections to rates for aligmnent with Brecht et al (2011) for NOx in the VTGCM
!  -- S. W. Bougher (170208)
!  Addition to permit NOVEM and NOINT to be calculated as VTGCM in Brecht et al (2011)
!  -- S. W. Bougher (170421)

  use ModSizeGitm
  use ModGITM
  use ModPlanet
  use ModRates
  use ModEUV
  use ModInputs, only: iDebugLevel, UseIonChemistry, UseNeutralChemistry,f107,f107a
  use ModConstants
  use ModSources, only: ChemicalHeatingS, IonPrecipIonRates,AuroralIonRates
  
  implicit None
    
    real :: Ions(nIons),Neutrals(nSpeciesTotal), eDensity
    real :: rtN2P_e, rtNOP_e, rtCO2P_e, rtO2P_e  ! DR rates
    real :: te, ti, tn
    real :: rtCO2_N2P, rtO_N2P, rtN2_OP, rtN4S_O, rtN4S_NO,  &
          rtN2D_O2                             ! Bi-molecular
  real :: rtO_O_CO2, rtO_CO_CO2, rtO_O2_CO2, rtO_N4S_CO2  ! Tri-molecular

! 11-Reactions and Rates taken from vtgcm  [May 2008] (T-independent)

  real :: rtCO2_OP, rtO_CO2P_tOP, rtO_CO2P_tO2P, rtN4S_O2P, rtN2D_CO2,  &
          rtN2D_CO, rtN2D_O, rtN2D_N2, rtN4S_CO2, rtNO_O2P, rtNO_N2D, &
          rtCO2P_O,rtCO2P_NO,rtCO2_N,rtCO2P_N2D,rtCOP_O,rtCOP_NO,&
          rtCOP_O2,rtCOP_CO2,rtO2P_N2D,rtO2P_N2D2,rtO2P_NO,rtO2P_N2,&
          rtN2P_N,rtN2P_CO,rtN2P_O2,rtN2P_O,rtN2P_NO,rtNP_CO2,&
          rtCO2_NP,rtOP_O2,rtOP_N2D,rtN4S_O2,rtN2D_NO,rtN2D_e,&
          rtNOP_e2,rtCO2P_N,rtCO2P_O2,rtO2P_N4S

  real :: tr05, tr104, expt1, expt2, rr_opn2
  real :: rt300te, rt300tn, rttn300
  logical :: useNeutralConstituent(nSpeciesTotal),useIonConstituent(nIons)
  integer :: nittot = 0
  real :: reactionrate(23)

contains
!-----------------------------------
subroutine calc_reaction_rates(iLon,iLat,iAlt,iBlock)
  
  integer, intent(in) :: iLon,iLat,iAlt,iBlock

            te = eTemperature(iLon,iLat,iAlt,iBlock)
           ti = iTemperature(iLon,iLat,iAlt,iBlock)
           tn = Temperature(iLon,iLat,iAlt,iBlock)*&
                TempUnit(iLon,iLat,iAlt)


!\ --------------------------------------------------------------------------
!  TN, TE, TI dependent reaction rates precalculated (like vtgcm)
!  Rates calculated once each iLON, iLAT, iALT bin of loop 
!  (m3/sec, m6/sec units needed) 
!  -- Only use Te , 1200 K rate for O2+ DR.  OK up to 200 km?
!/ --------------------------------------------------------------------------
! DR Rates with Electron Temp, te:  cm^3/s -> m^3/s
! real :: rtN2P_e, rtNOP_e, rtCO2P_e, rtO2P_e  ! DR rates

           rt300te = sqrt(300.0/te)
           rtN2P_e = 1.01e-7*(300./te)**0.39*1.0e-6
!          rtNOP_e = 3.4e-7*rt300te*1.0e-6
           rtNOP_e = 3.0e-7*rt300te*1.0e-6
!          rtNOP_e2 = 0.6e-07*rt300te*1.0e-6
           rtNOP_e2 = 1.0e-07*rt300te*1.0e-6
           rtCO2P_e = 3.5e-07*rt300te*1.0e-6
           rtO2P_e = 1.95e-07*(300./te)**0.70*1.0e-6
!\
! Reaction Rates with Neutral Temp, tn: cm^3/s -> m^3/s
! real :: rtCO2_N2P, rtO_N2P, rtN2_OP, rtN4S_O, rtN4S_NO, rtN2D_O2 ! Bi-molecular

           rt300tn = sqrt(300.0/tn)
           rttn300 = sqrt(tn/300.0)
           rtCO2_N2P = 9.0e-10*(300./ti)**0.23*1.0e-6
           rtO_N2P = 1.33e-10*(300./ti)**0.44*1.0e-6
           rtN2_OP = 1.20e-12*(300./ti)**0.45*1.0e-6
           rtN2D_O2 = 9.7e-12*exp(-185./tn)*1.0e-6
           rtN4S_O = 1.9e-17*rt300tn*(1.-0.57/tn**0.5)*1.e-06
           rtN4S_NO = 2.5e-10*rttn300*exp(-600./tn)*1.e-06
!\
! Reaction Rates with Neutral Temp, tn: cm^6/s -> m^6/s
! real :: rtO_O_CO2, rtO_CO_CO2, rtO_O2_CO2, rtO_N4S_CO2  ! Ter-molecular

!          rtO_O_CO2 = 2.75e-32*(200./tn)**3.3*1.e-12
           rtO_O_CO2 = 2.75e-32*1.e-12
           rtO_CO_CO2 = 6.5e-33*exp(-2180./tn)*1.e-12
!          rtO_O2_CO2 = 5.0e-28/tn**2.3*1.e-12
           rtO_O2_CO2 = 1.35e-33*1.e-12
!          rtO_N4S_CO2 = 2.0e-32*rt300tn*1.e-12
           rtO_N4S_CO2 = 1.83e-32*rt300tn*1.e-12
!\
! Reactions and Rates taken from vtgcm  [May 2008] (T-independent)
! real :: rtCO2_OP, rtO_CO2P_tOP, rtO_CO2P_tO2P, rtN4S_O2P, rtN2D_CO2,  &
!         rtN2D_CO, rtN2D_O, rtN2D_N2, rtN4S_CO2, rtNO_O2P, rtNO_N2D 
           
           rtCO2_OP = 1.1e-09*1.e-6 
           rtO_CO2P_tOP = 9.6e-11*1.e-6
           rtO_CO2P_tO2P = 1.64e-10*1.e-6 
           rtN4S_O2P = 1.0e-10*1.e-6 
!          rtN2D_CO2 = 3.6e-13*1.e-6 
           rtN2D_CO2 = 2.8e-13*1.e-6 
           rtN2D_CO = 1.9e-12*1.e-6 
!          rtN2D_O =  6.9e-13*1.e-6
           rtN2D_O =  2.0e-11*1.e-6
           rtN2D_N2 =  1.7e-14*1.e-6
           rtN4S_CO2 =  1.7e-16*1.e-6
           rtNO_O2P = 4.5e-10*1.e-6 
           rtNO_N2D =  6.7e-11*1.e-6

!\
!NEW RATES
!/

            rtCO2P_O2 = 5.5e-11*(300/Ti)**.82*1.0e-6
            rtCO2P_NO = 1.23e-10*1.e-6
            rtCO2P_N = 3.4e-10*1.e-6
            rtCO2P_N2D = 2.00e-10*1.0e-6
            rtCOP_O = 1.4e-10*1.0e-6
            rtCOP_NO = 4.2e-10*1.0e-6
            rtCOP_O2 = 1.5e-10*(300/Ti)**1.1*1.0e-6
            rtCOP_CO2 = 1.1e-9*1.e-6 
            rtO2P_N2D = 1.8e-10*1.e-6 
            rtO2P_N2D2 = 8.65e-11*1.e-6 
            rtO2P_NO = 4.5e-10*1.e-6 
            rtO2P_N2 = 1.0e-15*1.e-6 
            rtN2P_N = 1.0e-11*1.e-6 
            rtN2P_CO = 7.6e-11*1.e-6 
            rtN2P_O2 = 5.1e-11*(300/Ti)**1.16*1.e-6 
            rtN2P_O = 7.0e-12*(300/Ti)**.23*1.e-6 
            rtN2P_NO = 3.6e-10*1.e-6 
            rtNP_CO2 = 2.02e-10*1.e-6 
            rtCO2_NP = 9.18e-10*1.e-6 
            rtOP_O2 = 1.6e-11*(300/Ti)**.52*1.e-6 
            rtOP_N2D = 1.3e-10*1.e-6 
            rtN4S_O2 = 1.5e-14*Tn*exp(-3270/Tn)*1.e-6 
            rtN2D_NO = 6.7e-11*1.e-6 
            rtN2D_e = 3.86e-10*(Te/300)**.81*1.e-6 
            rtO2P_N4S = 1.0e-10*1.e-6

end subroutine calc_reaction_rates

subroutine calc_chemical_sources(iLon,iLat,iAlt,iBlock,IonSources, &
       IonLosses,NeutralSources, &
       NeutralLosses,ChemicalHeatingSub,Emission)


  integer, intent(in) :: ilon,ilat,ialt,iBlock
  real, intent(out) :: IonSources(nIons),IonLosses(nIons)
  real, intent(out) :: NeutralSources(nSpeciesTotal), NeutralLosses(nSpeciesTotal)
  real, intent(out) :: ChemicalHeatingSub,Emission(nEmissions)
  
  real :: Source, Reaction, rr
 
  integer :: iNeutral, iion

  
  real :: lon, lat, alt
  real, dimension(1:2) :: msis_temp
  real, dimension(1:8) :: msis_dens
  real :: LonDeg, LatDeg, AltKm, LST
  real,dimension(7)    :: AP  

  real :: percent,o2total,o2ptotal


   !---------------------------------------------------------------------------
 call start_timing("chemsources")
  if (.not.UseIonChemistry) return

  IonSources = 0.0
  NeutralSources = 0.0
  IonLosses  = 0.0
  NeutralLosses = 0.0
  
  ChemicalHeatingSub = 0.0
  Emission = 0.0

  o2ptotal = 0
  o2total = 0


!\
! Nitrogen Photochemistry:--------------------------------------------------+
!/
              ! ----------------------------------------------------------
              ! N2 + hv ==> N(4S) + N(2D)   Assume 50% Branching Ratio
              ! N2 + pe ==> N(4S) + N(2D)   Assume 50% Branching Ratio (later)
              ! Triple photodissrate due to fine structure approximation (VTGCM)
              ! ----------------------------------------------------------
              !rr=EuvDissRateS(iLon,iLat,iAlt,iN2_,iBlock) + PEDissRateS(iLon,iLat,iAlt,iN2_,iBlock)
              !rr=EuvDissRateS(iLon,iLat,iAlt,iN2_,iBlock)
              rr=EuvDissRateS(iLon,iLat,iAlt,iN2_,iBlock)*3.0 

              Reaction = rr * Neutrals(iN2_)

              NeutralLosses(iN2_) = NeutralLosses(iN2_) + Reaction
              NeutralSources(iN4S_) = NeutralSources(iN4S_) + 2.*0.5*Reaction
              NeutralSources(iN2D_) = NeutralSources(iN2D_) + 2.*0.5*Reaction

              ! ----------------------------------------------------------
              ! N2 + hv ==> N2+ 
              ! ----------------------------------------------------------

              Reaction = EuvIonRateS(iLon,iLat,iAlt,iN2P_,iBlock)*Neutrals(iN2_)

              NeutralLosses(iN2_) = NeutralLosses(iN2_) + Reaction
              IonSources(iN2P_) = IonSources(iN2P_) + Reaction
              reactionrate(10) = reaction
!\
! CO2 Photochemistry:-----------------------------------------------------+
!/
              ! ----------------------------------------------------------
              ! CO2 + hv ==> CO2+ 
              ! ----------------------------------------------------------

              Reaction = EuvIonRateS(iLon,iLat,iAlt,iCO2P_,iBlock)*Neutrals(iCO2_)
!if (ialt .eq. 42) then
!   write(*,*) reaction,EuvIonRateS(iLon,iLat,iAlt,iCO2P_,iBlock),Neutrals(iCO2_)
!endif
             NeutralLosses(iCO2_) = NeutralLosses(iCO2_) + Reaction
!              NeutralLosses(iCO2_) = 0.0
              IonSources(iCO2P_) = IonSources(iCO2P_) + Reaction

!              write(*,*) "in calc chem: ",reaction,euvionrates(1,1,ialt,ico2P_,iblock), neutrals(ico2_)
              reactionrate(1) = reaction
              !if (ialt .eq. 50) write(*,*) " CO2 + hv ==> CO2+  ",niters,&
!                   ialt,reaction,euvionrates(1,1,ialt,ico2p_,iblock),neutrals(ico2_)
!write(*,*) reaction,EuvIonRateS(iLon,iLat,iAlt,iCO2P_,iBlock),Neutrals(iCO2_)
!stop
              ! ----------------------------------------------------------
              ! CO2 + hv ==> CO + O 
              ! ----------------------------------------------------------

              Reaction = EuvDissRateS(iLon,iLat,iAlt,iCO2_,iBlock)*Neutrals(iCO2_)

             NeutralLosses(iCO2_) = NeutralLosses(iCO2_) + Reaction
!              NeutralLosses(iCO2_) = 0.0
              NeutralSources(iO_) = NeutralSources(iO_) + Reaction
              NeutralSources(iCO_) = NeutralSources(iCO_) + Reaction
!\
! O2 Photochemistry:-----------------------------------------------------+
!
!             ! ----------------------------------------------------------
!             ! O2 + hv ==> O2+   Minor source of O2+  (Mnor)
!             ! ----------------------------------------------------------
!
              Reaction = EuvIonRateS(iLon,iLat,iAlt,iO2P_,iBlock)*Neutrals(iO2_)
 
              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction
              IonSources(iO2P_) = IonSources(iO2P_) + Reaction
              !/
              reactionrate(11) = reaction

              ! ----------------------------------------------------------
              ! O2 + hv ==> O + O 
              ! ----------------------------------------------------------

              Reaction = EuvDissRateS(iLon,iLat,iAlt,iO2_,iBlock)*Neutrals(iO2_)

              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction
              NeutralSources(iO_) = NeutralSources(iO_) + 2.*Reaction
!\
! O Photochemistry:-----------------------------------------------------+
!/
              ! ----------------------------------------------------------
              ! O + hv ==> O+ 
              !  Two Components:
              !  (1)  O   --> O+ direct (total ionization): High Altitude
              !       photoion(1:NWH,iOP_) slot filled in fill_photo
              !  (2)  CO2 --> O+ indirect (branching to O+ ionization): Near F1-peak
              !       photoion(1:NWH,iNOP_) slot filled in fill_photo
              ! ----------------------------------------------------------

              Reaction = EuvIonRateS(iLon,iLat,iAlt,iOP_,iBlock)*Neutrals(iO_) + &
                         EuvIonRateS(iLon,iLat,iAlt,iNOP_,iBlock)*Neutrals(iCO2_)

              NeutralLosses(iO_) = NeutralLosses(iO_) + Reaction
              IonSources(iOP_) = IonSources(iOP_) + Reaction
              reactionrate(6) = reaction
!              if (ialt .eq. 41) write(*,*) "O+hv: ",reaction
!\
!-----------End Photochemistry-----------------------------------------------------+
!\
!-----------Begin Electron Recombination Chemistry---------------------------------+
!/
              ! -----------------------------------------------------------
              ! N2+ + e- ==> 2.*(0.25*N4S + 0.75*N2D)  Fix and Sung (2001) Branching
              ! -----------------------------------------------------------

              Reaction = rtN2P_e * Ions(ie_) * Ions(iN2P_) 

              IonLosses(iN2P_) = IonLosses(iN2P_) + Reaction
              NeutralSources(iN4S_) = NeutralSources(iN4S_) + 2.*0.25*Reaction
              NeutralSources(iN2D_) = NeutralSources(iN2D_) + 2.*0.75*Reaction
              reactionrate(12) = reaction
              ! -----------------------------------------------------------
              ! O2+ + e- ==> 2O
              ! -----------------------------------------------------------

              Reaction = rtO2P_e * Ions(ie_) * Ions(iO2P_) 

              IonLosses(iO2P_) = IonLosses(iO2P_) + Reaction
              NeutralSources(iO_) = NeutralSources(iO_) + 2.*Reaction
              reactionrate(13) = reaction
              ! -----------------------------------------------------------
              ! CO2+ + e- ==> O + CO
              ! -----------------------------------------------------------

              Reaction =rtCO2P_e * Ions(ie_) * Ions(iCO2P_) 
              reactionrate(2) = reaction
              IonLosses(iCO2P_) = IonLosses(iCO2P_) + Reaction
              NeutralSources(iO_) = NeutralSources(iO_) + Reaction
              NeutralSources(iCO_) = NeutralSources(iCO_) + Reaction
              !if (ialt .eq. 50) write(*,*) "CO2+ + e- ==> O + CO",niters, ialt,reaction,rtCO2P_e , Ions(ie_) , Ions(iCO2P_) 
               ! -----------------------------------------------------------
              ! NO+ + e- ==> O + N2D 
              ! g = 0.75
              ! -----------------------------------------------------------

              Reaction = rtNOP_e * Ions(ie_) * Ions(iNOP_) 
              !if (ialt .eq. 47)               write(*,*) -1*reaction, "NO+ + e- ==> O + N2D ",rtNOP_e,Ions(ie_) ,Ions(iNOP_) 
              IonLosses(iNOP_) = IonLosses(iNOP_) + Reaction
              NeutralSources(iN2D_) = NeutralSources(iN2D_) + Reaction
              NeutralSources(iO_) = NeutralSources(iO_) + Reaction
              reactionrate(14) = reaction

              ! -----------------------------------------------------------
              ! NO+ + e- ==> O + N4S 
              ! (1-g) = 0.25
              ! -----------------------------------------------------------

              Reaction = rtNOP_e2 * Ions(ie_) * Ions(iNOP_) 
              !if (ialt .eq. 47)  write(*,*) -1*reaction, "NO+ + e- ==> O + N4S ",rtNOP_e2,Ions(ie_) ,Ions(iNOP_) 
              IonLosses(iNOP_) = IonLosses(iNOP_) + Reaction
              NeutralSources(iN4S_) = NeutralSources(iN4S_) + Reaction
              NeutralSources(iO_) = NeutralSources(iO_) + Reaction

!-----------End Electron Recombination Chemistry-------------------------------------------+
!/
!-----------Begin Bi-Molecular Ion-Neutral and Neutral-Neutral Chemistry-----------+
!/
              ! -----------------------------------------------------------
              ! CO2+ + O ==> O2+ + CO  Fast 
              ! -----------------------------------------------------------
              Reaction = rtO_CO2P_tO2P * Neutrals(iO_) * Ions(iCO2P_) 
              reactionrate(3) = reaction
!if (ialt .eq. 42) then
!   write(*,*) reaction,rtO_CO2P_tO2P,Neutrals(iO_) ,Ions(iCO2P_) 


!endif
              NeutralLosses(iO_) = NeutralLosses(iO_) + Reaction
              IonLosses(iCO2P_) = IonLosses(iCO2P_) + Reaction
              NeutralSources(iCO_) = NeutralSources(iCO_) + Reaction
              IonSources(iO2P_) = IonSources(iO2P_) + Reaction

              !if (ialt .eq. 50) write(*,*) "CO2+ + O ==> O2+ + CO",niters,&
!                   ialt,reaction,rtO_CO2P_tO2P,  Neutrals(iO_) , Ions(iCO2P_) 
              ! -----------------------------------------------------------
              ! CO2+ + O ==> O+ + CO2  Fast 
              ! -----------------------------------------------------------
              Reaction =rtO_CO2P_tOP * Neutrals(iO_) * Ions(iCO2P_) 
!if (ialt .eq. 42) then
!   write(*,*) reaction,rtO_CO2P_tOP,Neutrals(iO_) ,Ions(iCO2P_) 
!stop
!endif
              reactionrate(4) = reaction
              NeutralLosses(iO_) = NeutralLosses(iO_) + Reaction
              IonLosses(iCO2P_) = IonLosses(iCO2P_) + Reaction
             NeutralSources(iCO2_) = NeutralSources(iCO2_) + Reaction
!              NeutralSources(iCO2_) = 0.0
              IonSources(iOP_) = IonSources(iOP_) + Reaction
              reactionrate(7) = reaction
              !if (ialt .eq. 50) write(*,*) "CO2+ + O ==> O+ + CO2",niters, &
              !ialt,reaction ,rtO_CO2P_tOP ,Neutrals(iO_) , Ions(iCO2P_) 
!              if (ialt .eq. 41)              write(*,*) "CO2+ + O ",reaction, neutrals(iO_)
              ! -----------------------------------------------------------
              ! CO2 + O+ ==> O2+ + CO  Fast 
              ! -----------------------------------------------------------
              Reaction = rtCO2_OP * Neutrals(iCO2_) * Ions(iOP_) 

!                write(*,*)ialt, niters, rtCO2_OP ,neutrals(iCO2_),ions(iOP_)
              !if (ialt .eq. 45 .or. ialt  .eq. 46)
!              write(*,*)"Before: ", ialt,rtco2_op,Neutrals(iCO2_),&
!                   Ions(iOP_) ,reaction
             NeutralLosses(iCO2_) = NeutralLosses(iCO2_) + Reaction
!              NeutralLosses(iCO2_) = 0.0
              IonLosses(iOP_) = IonLosses(iOP_) + Reaction
              NeutralSources(iCO_) = NeutralSources(iCO_) + Reaction
              IonSources(iO2P_) = IonSources(iO2P_) + Reaction
              reactionrate(8) = reaction
!              write(*,*)"After: ", ialt,rtco2_op,Neutrals(iCO2_),&
!                   Ions(iOP_)-reaction ,reaction
!              if (ialt .eq. 41) write(*,*) "CO2 + O+: ",reaction, neutrals(ico2_)
              ! -----------------------------------------------------------
              ! CO2 + N2+ ==> CO2+ + N2  
              ! -----------------------------------------------------------
              Reaction = rtCO2_N2P * Neutrals(iCO2_) * Ions(iN2P_) 
              reactionrate(5) = reaction
              NeutralLosses(iCO2_) = NeutralLosses(iCO2_) + Reaction
!              NeutralLosses(iCO2_) = 0.0
              IonLosses(iN2P_) = IonLosses(iN2P_) + Reaction
              NeutralSources(iN2_) = NeutralSources(iN2_) + Reaction
              IonSources(iCO2P_) = IonSources(iCO2P_) + Reaction

                           
              ! -----------------------------------------------------------
              ! O + N2+ ==> NO+ + N4S  
              ! -----------------------------------------------------------
              Reaction = rtO_N2P * Neutrals(iO_) * Ions(iN2P_) 
              !if (ialt .eq. 47) write(*,*) reaction, "O + N2+ ==> NO+ + N4S"
              NeutralLosses(iO_) = NeutralLosses(iO_) + Reaction
              IonLosses(iN2P_) = IonLosses(iN2P_) + Reaction
              NeutralSources(iN4S_) = NeutralSources(iN4S_) + Reaction
              IonSources(iNOP_) = IonSources(iNOP_) + Reaction
              reactionrate(16) = reaction
              ! -----------------------------------------------------------
              ! N2 + O+ ==> NO+ + N4S  
              ! -----------------------------------------------------------
              Reaction = rtN2_OP * Neutrals(iN2_) * Ions(iOP_) 
              !if (ialt .eq. 47) write(*,*) reaction, "N2 + O+ ==> NO+ + N4S "
              NeutralLosses(iN2_) = NeutralLosses(iN2_) + Reaction
              IonLosses(iOP_) = IonLosses(iOP_) + Reaction
              NeutralSources(iN4S_) = NeutralSources(iN4S_) + Reaction
              IonSources(iNOP_) = IonSources(iNOP_) + Reaction
              reactionrate(17) = reaction
!              if (ialt .eq. 41) write(*,*) "N2 + O+: ",reaction
!              if (ialt .eq. 41) stop
              ! -----------------------------------------------------------
              ! N4S+ O2+ ==> NO+ + O  
              ! -----------------------------------------------------------
              Reaction = rtN4S_O2P * Neutrals(iN4S_) * Ions(iO2P_) 
              !if (ialt .eq. 47) write(*,*) reaction, "N4S+ O2+ ==> NO+ + O  "
              NeutralLosses(iN4S_) = NeutralLosses(iN4S_) + Reaction
              IonLosses(iO2P_) = IonLosses(iO2P_) + Reaction
              NeutralSources(iO_) = NeutralSources(iO_) + Reaction
              IonSources(iNOP_) = IonSources(iNOP_) + Reaction
              reactionrate(18) = reaction

!----------------------------------
!New reactions
! ----------------------------------
!               ! -----------------------------------------------------------
!               ! CO2+ + O2 ==> CO2 + O2+
!               ! -----------------------------------------------------------
!               reaction = rtCO2P_O2 * Neutrals(iO2_) * Ions(iCO2P_)
!               NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction
!               IonLosses(iCO2P_) = IonLosses(iCO2P_) + Reaction
!               NeutralSources(iCO2_) = NeutralSources(iCO2_) + Reaction
!               IonSources(iO2P_) = IonSources(iO2P_) + Reaction
!                
!                
!                ! -----------------------------------------------------------
!                !CO2+ + NO ==> NO+ + CO2
!                ! -----------------------------------------------------------
!                reaction = rtCO2P_NO * Neutrals(iNO_) * Ions(iCO2P_)
!               !if (ialt .eq. 47)  write(*,*) reaction, "CO2+ + NO ==> NO+ + CO2"
!                NeutralLosses(iNO_) = NeutralLosses(iNO_) + Reaction
!                IonLosses(iCO2P_) = IonLosses(iCO2P_) + Reaction
!                NeutralSources(iCO2_) = NeutralSources(iCO2_) + Reaction
!                IonSources(iNOP_) = IonSources(iNOP_) + Reaction
!                reactionrate(19) = reaction
!  
!                 ! -----------------------------------------------------------
!                ! CO2+ + N2D ==> N+ + CO2
!                ! -----------------------------------------------------------
!  !              reaction = rtCO2P_N2S * Neutrals(iN2D_) * Ions(iCO2P_)
!  !              NeutralLosses(iN2D_) = NeutralLosses(iN2D_) + Reaction
!  !              IonLosses(iCO2P_) = IonLosses(iCO2P_) + Reaction
!  !              NeutralSources(iCO2_) = NeutralSources(iCO2_) + Reaction
!  !              IonSources(iNP_) = IonSources(iNP_) + Reaction
!  
!  
!                ! -----------------------------------------------------------
!                ! O2+ + N4S ==> NO+ + O
!                ! -----------------------------------------------------------
!                reaction = rtO2P_N4S * Neutrals(iN4S_) * Ions(iO2P_)
!                !if (ialt .eq. 47) write(*,*) reaction, "O2+ + N4S ==> NO+ + O"
!                NeutralLosses(iN4S_) = NeutralLosses(iN4S_) + Reaction
!                IonLosses(iO2P_) = IonLosses(iO2P_) + Reaction
!                NeutralSources(iO_) = NeutralSources(iO_) + Reaction
!                IonSources(iNOP_) = IonSources(iNOP_) + Reaction
!                reactionrate(20) = reaction
!                ! -----------------------------------------------------------
!                ! O2+ + N2D ==> NO+ + O
!                ! -----------------------------------------------------------
!                reaction = rtO2P_N2D * Neutrals(iN2D_) * Ions(iO2P_)
!                !if (ialt .eq. 47) write(*,*) reaction, "O2+ + N2D ==> NO+ + O"
!                NeutralLosses(iN2D_) = NeutralLosses(iN2D_) + Reaction
!                IonLosses(iO2P_) = IonLosses(iO2P_) + Reaction
!                NeutralSources(iO_) = NeutralSources(iO_) + Reaction
!                IonSources(iNOP_) = IonSources(iNOP_) + Reaction
!                reactionrate(21) = reaction
!                ! -----------------------------------------------------------
!                ! O2+ + N2 ==> NO+ + NO
!                ! ----------------------------------------------------------- 
!                reaction = rtO2P_N2 * Neutrals(iN2_) * Ions(iO2P_)
!               !if (ialt .eq. 47)  write(*,*) reaction, "O2+ + N2 ==> NO+ + NO",rtO2P_N2,Neutrals(iN2_) ,Ions(iO2P_)
!                NeutralLosses(iN2_) = NeutralLosses(iN2_) + Reaction
!                IonLosses(iO2P_) = IonLosses(iO2P_) + Reaction
!                NeutralSources(iNO_) = NeutralSources(iNO_) + Reaction
!                IonSources(iNOP_) = IonSources(iNOP_) + Reaction
!                reactionrate(22) = reaction
!  
!  
!                ! -----------------------------------------------------------
!                ! N2+ + O2 ==> N2 + O2+
!                ! ----------------------------------------------------------- 
!                reaction = rtN2P_O2 * Neutrals(iO2_) * Ions(iN2P_)
!                NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction
!                IonLosses(iN2P_) = IonLosses(iN2P_) + Reaction
!                NeutralSources(iN2_) = NeutralSources(iN2_) + Reaction
!                IonSources(iO2P_) = IonSources(iO2P_) + Reaction
!  
!                ! -----------------------------------------------------------
!                ! N2+ + O ==> O+ + N2
!                ! ----------------------------------------------------------- 
!                reaction = rtN2P_O * Neutrals(iO_) * Ions(iN2P_)
!                NeutralLosses(iO_) = NeutralLosses(iO_) + Reaction
!                IonLosses(iN2P_) = IonLosses(iN2P_) + Reaction
!                NeutralSources(iN2_) = NeutralSources(iN2_) + Reaction
!                IonSources(iOP_) = IonSources(iOP_) + Reaction
!  
!                ! -----------------------------------------------------------
!                ! N2+ + NO ==> N2 + NO+
!                ! ----------------------------------------------------------- 
!                reaction = rtN2P_NO* Neutrals(iNO_) * Ions(iN2P_)
!                !if (ialt .eq. 47) write(*,*) reaction, "N2+ + NO ==> N2 + NO+"
!                NeutralLosses(iNO_) = NeutralLosses(iNO_) + Reaction
!                IonLosses(iN2P_) = IonLosses(iN2P_) + Reaction
!                NeutralSources(iN2_) = NeutralSources(iN2_) + Reaction
!                IonSources(iNOP_) = IonSources(iNOP_) + Reaction
!                reactionrate(23) = reaction
!  
!                 ! -----------------------------------------------------------
!                ! O+ + O2 ==> O + O2+
!                ! ----------------------------------------------------------- 
!                reaction = rtOP_O2* Neutrals(iO2_) * Ions(iOP_)
!                NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction
!                IonLosses(iOP_) = IonLosses(iOP_) + Reaction
!                NeutralSources(iO_) = NeutralSources(iO_) + Reaction
!                IonSources(iO2P_) = IonSources(iO2P_) + Reaction

!\
!-----------End Bi-Molecular Ion-Neutral and Neutral-NeutralChemistry-------------------+
!/
!\
!-----------Ter-Molecular Neutral-Neutral Chemistry-------------------------------+
!/
              ! -----------------------------------------------------------
              ! O + O + CO2 ==> O2 + CO2  
              ! -----------------------------------------------------------
              Reaction = rtO_O_CO2*Neutrals(iO_)*Neutrals(iO_)*Neutrals(iCO2_)

              NeutralLosses(iO_) = NeutralLosses(iO_) + 2.*Reaction
              NeutralSources(iO2_) = NeutralSources(iO2_) + Reaction

              ! -----------------------------------------------------------
              ! O + O2 + CO2 ==> stuff + CO2  
              ! -----------------------------------------------------------
!             Reaction = rtO_O_CO2 *Neutrals(iO_)*Neutrals(iO_)*Neutrals(iCO2_)
              Reaction = rtO_O2_CO2 *Neutrals(iO_)*Neutrals(iO2_)*Neutrals(iCO2_)

              NeutralLosses(iO_) = NeutralLosses(iO_) + Reaction
              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction

              ! -----------------------------------------------------------
              ! O + CO + CO2 ==> 2.*CO2  
              ! -----------------------------------------------------------
              Reaction = rtO_CO_CO2 *Neutrals(iO_)*Neutrals(iCO_)*Neutrals(iCO2_)

              NeutralLosses(iO_) = NeutralLosses(iO_) + Reaction
              NeutralLosses(iCO_) = NeutralLosses(iCO_) + Reaction
             NeutralSources(iCO2_) = NeutralSources(iCO2_) + Reaction
!              NeutralSources(iCO2_) = 0.0

              ! -----------------------------------------------------------
              ! O + N4S + CO2 ==> NO + CO2
              ! -----------------------------------------------------------
              Reaction = rtO_N4S_CO2*Neutrals(iO_)*Neutrals(iN4S_)*Neutrals(iCO2_)

              NeutralLosses(iO_) = NeutralLosses(iO_) + Reaction
              NeutralLosses(iN4S_) = NeutralLosses(iN4S_) + Reaction
              NeutralSources(iNO_) = NeutralSources(iNO_) + Reaction
!\
!-----------NOX Specific Ion-Neutral and Neutral-Neutral Chemistry-------------------+
!/
              ! -----------------------------------------------------------
              ! N2D + CO2  ==> NO + CO
              ! -----------------------------------------------------------
              Reaction = rtN2D_CO2*Neutrals(iN2D_)*Neutrals(iCO2_)

              NeutralLosses(iN2D_) = NeutralLosses(iN2D_) + Reaction
             NeutralLosses(iCO2_) = NeutralLosses(iCO2_) + Reaction
!              NeutralLosses(iCO2_) = 0.0
              NeutralSources(iNO_) = NeutralSources(iNO_) + Reaction
              NeutralSources(iCO_) = NeutralSources(iCO_) + Reaction
 
              ! -----------------------------------------------------------
              ! N2D + CO  ==> N4S + CO
              ! -----------------------------------------------------------
              Reaction = rtN2D_CO*Neutrals(iN2D_)*Neutrals(iCO_)

              NeutralLosses(iN2D_) = NeutralLosses(iN2D_) + Reaction
              NeutralSources(iN4S_) = NeutralSources(iN4S_) + Reaction

              ! -----------------------------------------------------------
              ! N2D + O  ==> N4S + O
              ! -----------------------------------------------------------
              Reaction = rtN2D_O*Neutrals(iN2D_)*Neutrals(iO_)

              NeutralLosses(iN2D_) = NeutralLosses(iN2D_) + Reaction
              NeutralSources(iN4S_) = NeutralSources(iN4S_) + Reaction

              ! -----------------------------------------------------------
              ! N2D + O2  ==> N4S + O2
              ! -----------------------------------------------------------
              Reaction = rtN2D_O2*Neutrals(iN2D_)*Neutrals(iO2_)

              NeutralLosses(iN2D_) = NeutralLosses(iN2D_) + Reaction
              NeutralSources(iN4S_) = NeutralSources(iN4S_) + Reaction

              ! -----------------------------------------------------------
              ! N2D + N2  ==> N4S + N2
              ! -----------------------------------------------------------
              Reaction = rtN2D_N2*Neutrals(iN2D_)*Neutrals(iN2_)

              NeutralLosses(iN2D_) = NeutralLosses(iN2D_) + Reaction
              NeutralSources(iN4S_) = NeutralSources(iN4S_) + Reaction

              ! -----------------------------------------------------------
              ! NO + N4S  ==> N2 + O
              ! -----------------------------------------------------------
              Reaction = rtN4S_NO*Neutrals(iN4S_)*Neutrals(iNO_)

              NeutralLosses(iNO_) = NeutralLosses(iNO_) + Reaction
              NeutralLosses(iN4S_) = NeutralLosses(iN4S_) + Reaction
              NeutralSources(iN2_) = NeutralSources(iN2_) + Reaction
              NeutralSources(iO_) = NeutralSources(iO_) + Reaction

              ! -----------------------------------------------------------
              ! N4S + O  ==> NO* + hv (Nightglow reaction)
              ! New code added by S. W. Bougher (170427)
              ! -----------------------------------------------------------
              Reaction = rtN4S_O*Neutrals(iN4S_)*Neutrals(iO_)

              NeutralLosses(iN4S_) = NeutralLosses(iN4S_) + Reaction
              NeutralLosses(iO_) = NeutralLosses(iO_) + Reaction
              NeutralSources(iNO_) = NeutralSources(iNO_) + Reaction

!             Total NO (delta plus gamma band system) mission (190-260 nm)
!             Initially Should be ph/m3.sec units  (check)
!             Convert to ph/cm3.sec units  
              Emission(iENOUV_) = Emission(iENOUV_) + Reaction*1.0E-06

!----------------------------------
!New reactions
!----------------------------------
!               ! -----------------------------------------------------------
!               ! N + CO2  ==> NO + CO
!               ! -----------------------------------------------------------
!               Reaction = rtN4S_CO2*Neutrals(iN4S_)*Neutrals(iCO2_)
! 
!               NeutralLosses(iN4S_) = NeutralLosses(iN4S_) + Reaction
!               NeutralLosses(iCO2_) = NeutralLosses(iCO2_) + Reaction
!               NeutralSources(iNO_) = NeutralSources(iNO_) + Reaction
!               NeutralSources(iCO_) = NeutralSources(iCO_) + Reaction
! 
!               ! -----------------------------------------------------------
!               ! N(2D) + NO  ==> N2 + O
!               ! -----------------------------------------------------------
!               Reaction = rtN2D_NO*Neutrals(iN2D_)*Neutrals(iNO_)
! 
!               NeutralLosses(iN2D_) = NeutralLosses(iN2D_) + Reaction
!               NeutralLosses(iNO_) = NeutralLosses(iNO_) + Reaction
!               NeutralSources(iN2_) = NeutralSources(iN2_) + Reaction
!               NeutralSources(iO_) = NeutralSources(iO_) + Reaction  
! 
!               ! -----------------------------------------------------------
!               ! N(2D) + e  ==> N(4S) + e
!               ! -----------------------------------------------------------
!               Reaction = rtN2D_e*Neutrals(iN2D_)*Ions(ie_)
! 
!               NeutralLosses(iN2D_) = NeutralLosses(iN2D_) + Reaction
!               NeutralSources(iN4S_) = NeutralSources(iN4S_) + Reaction
!\
!-----------Inert Species (zeroed out sources and sinks)
!-----------Species not calculated yet (zeroed out sources and sinks)
!/
              ! -----------------------------------------------------------
              ! He, H:  little idea about Aij coefficients
              ! -----------------------------------------------------------
!------------CO2 infinite reservoir assumption used-------------
              NeutralLosses(iCO2_) = 0.0
              NeutralSources(iCO2_) = 0.0
!------------Diffusion and advection only (no chemistry presently) -------
              NeutralLosses(iAr_) = 0.0
              NeutralSources(iAr_) = 0.0
              NeutralLosses(iHe_) = 0.0
              NeutralSources(iHe_) = 0.0
              NeutralLosses(iH_) = 0.0
              NeutralSources(iH_) = 0.0
!\
!------------------------------------------------------+
!-----------End Total Chemistry-------------------------------------------+
!------------------------------------------------------+

if (.NOT. useNeutralChemistry) then 

   NeutralSources = 0.0
   NeutralLosses = 0.0
endif

call end_timing("chemsources")    
   end subroutine calc_chemical_sources



end Module ModChemistry

