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
  real, dimension(1,3) :: Vibrational_Population
  real, dimension(3,4) :: Branching_Per_Vibration = &
       transpose(reshape((/ 0.265, 0.473, 0.204, 0.058, &
                            0.073, 0.278, 0.51, 0.139,  &
                            0.02, 0.764, 0.025, 0.211 /), &
                            (/4,3/)))
  real, dimension(4,1) :: Exotherm_O2_DR
  real, dimension(1,1) :: avgHeatMatrix
  real :: avgHeat
  

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
!  -- Only use Te <= 1200 K rate for O2+ DR.  OK up to ~220 km
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
           rtO2P_e = 2.4e-07*(300./te)**0.70*1.0e-6
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
           real :: beta = 0.75
           real :: quenchingFactor = 0.0
           real :: fO2IR = 0.0
           real :: o2ver = 0.0
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
            !!rr=EuvDissRateS(iLon,iLat,iAlt,iN2_,iBlock) + PEDissRateS(iLon,iLat,iAlt,iN2_,iBlock)
            !!rr=EuvDissRateS(iLon,iLat,iAlt,iN2_,iBlock)
            rr=EuvDissRateS(iLon,iLat,iAlt,iN2_,iBlock)*3.0 

            Reaction = rr
 
            NeutralLosses(iN2_) = NeutralLosses(iN2_) + Reaction
            NeutralSources(iN4S_) = NeutralSources(iN4S_) + 2.*0.5*Reaction
            NeutralSources(iN2D_) = NeutralSources(iN2D_) + 2.*0.5*Reaction

            !No chemical heating included in this reaction for Earth GITM
            
            ! ----------------------------------------------------------
            ! N2 + hv ==> N2+ 
            ! ----------------------------------------------------------
            !BP: Currently there are no loss terms for the N2 ionization so it is removed

            !Reaction = EuvIonRateS(iLon,iLat,iAlt,iN2P_,iBlock)

            !NeutralLosses(iN2_) = NeutralLosses(iN2_) + Reaction
            !IonSources(iN2P_) = IonSources(iN2P_) + Reaction
            !reactionrate(10) = reaction

           !\
           ! CO2 Photochemistry:-----------------------------------------------------+
           !/
           ! ----------------------------------------------------------
           ! CO2 + hv ==> CO2+
           ! ----------------------------------------------------------

           Reaction = EuvIonRateS(iLon,iLat,iAlt,iCO2P_,iBlock)
          
           NeutralLosses(iCO2_) = NeutralLosses(iCO2_) + Reaction
           IonSources(iCO2P_) = IonSources(iCO2P_) + Reaction

           !reactionrate(1) = reaction
           !write(*,*) reactionrate(1)

           ! ----------------------------------------------------------
           ! CO2 + hv (5.453 eV) ==> CO + O
           ! Dissociation threshold from Shaw et al., 1995
           ! ----------------------------------------------------------

           Reaction = EuvDissRateS(iLon,iLat,iAlt,iCO2_,iBlock)

           NeutralLosses(iCO2_) = NeutralLosses(iCO2_) + Reaction
           !!   NeutralLosses(iCO2_) = 0.0
           NeutralSources(iO_) = NeutralSources(iO_) + Reaction
           NeutralSources(iCO_) = NeutralSources(iCO_) + Reaction
          
           !reactionrate(2) = reaction

           ! ----------------------------------------------------------         
           ! O + hv ==> O+                  
           ! ----------------------------------------------------------

           Reaction = EuvIonRateS(iLon,iLat,iAlt,iOP_,iBlock)

           NeutralLosses(iO_) = NeutralLosses(iO_) + Reaction
           IonSources(iOP_) = IonSources(iOP_) + Reaction
           !reactionrate(3) = reaction

           !No chemical heating included in the Earth GITM code for
           !this ionization

           ! ----------------------------------------------------------
           ! CO2+ + O ==> O2+ + CO + 1.33 eV
           ! According to Fox and Hac, 2009
           ! ---------------------------------------------------------- 

           Reaction = rtO_CO2P_tO2P * Neutrals(iO_) * Ions(iCO2P_) 
           !reactionrate(4) = reaction

           NeutralLosses(iO_) = NeutralLosses(iO_) + Reaction
           IonLosses(iCO2P_) = IonLosses(iCO2P_) + Reaction
           NeutralSources(iCO_) = NeutralSources(iCO_) + Reaction
           IonSources(iO2P_) = IonSources(iO2P_) + Reaction

           ChemicalHeatingSub = ChemicalHeatingSub + &
                Reaction * 1.33

           ! -----------------------------------------------------------
           ! CO2+ + O ==> O+ + CO2  Fast 
           ! -----------------------------------------------------------
           Reaction =rtO_CO2P_tOP * Neutrals(iO_) * Ions(iCO2P_) 

           NeutralLosses(iO_) = NeutralLosses(iO_) + Reaction
           IonLosses(iCO2P_) = IonLosses(iCO2P_) + Reaction
           NeutralSources(iCO2_) = NeutralSources(iCO2_) + Reaction
           IonSources(iOP_) = IonSources(iOP_) + Reaction
           
           !reactionrate(5) = reaction

           ! -----------------------------------------------------------
           ! CO2 + O+ ==> O2+ + CO + 1.21 eV Fast
           ! See Gu et al., 2020 paper for 1.21 eV, but may be 0.65 eV to KE
           ! Need to confirm with Aaron/Steve
           ! -----------------------------------------------------------
           Reaction = rtCO2_OP * Neutrals(iCO2_) * Ions(iOP_) 

           NeutralLosses(iCO2_) = NeutralLosses(iCO2_) + Reaction

           IonLosses(iOP_) = IonLosses(iOP_) + Reaction
           NeutralSources(iCO_) = NeutralSources(iCO_) + Reaction
           IonSources(iO2P_) = IonSources(iO2P_) + Reaction

           ChemicalHeatingSub = ChemicalHeatingSub + &
                Reaction * 1.21
           ! -----------------------------------------------------------
           ! O2+ + e- ==> 2O
           ! -----------
              ! O2+ + e -> O(3P) + O(3P) + 6.99 eV
              !         -> O(3P) + O(1D) + 5.02 eV
              !         -> O(1D) + O(1D) + 3.06 eV
              !         -> O(1D) + O(1S) + 0.83 eV
           ! -----------
           !This needs work. I should do a varying system
           !for the vibrational states
           !                  nu = 0  = (0.265, 0.473, 0.204, 0.058)
           !                  nu = 1  = (0.073, 0.278, 0.51, 0.139)
           !                  nu = 2  = (0.02, 0.764, 0.025, 0.211)
           ! -----------------------------------------------------------

           Reaction = rtO2P_e * Ions(ie_) * Ions(iO2P_) 

           IonLosses(iO2P_) = IonLosses(iO2P_) + Reaction
           NeutralSources(iO_) = NeutralSources(iO_) + 2.*Reaction

           Exotherm_O2_DR(:,1) = (/6.99, 5.02, 3.06, 0.83/) !eV

           !Find the correct fractional population of the vibration states
           !based on altitude
           if (Altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0 .le. 130.0) then
              !All in the nu = 0 state (see Fig. 1 of Fox 1985, THE 02+ VIB-
              !RATIONAL DISTRIBUTION IN THE VENUSIAN IONOSPHERE)
             Vibrational_Population(1,:) = (/1, 0, 0/)
           else
              !Linearly decreasing system from 130km - 170 km
              !just y = mx + b from Figure 1 for nu = 0, 1 and 2
              Vibrational_Population(1,1) = &
                   -0.27/(40*10**3) * (Altitude_GB(iLon,iLat,iAlt,iBlock) - 130*10**3) + 1.0
              Vibrational_Population(1,2) = &
                   0.07/(40*10**3) * (Altitude_GB(iLon,iLat,iAlt,iBlock) - 130*10**3)
              Vibrational_Population(1,3) = &
                   0.06/(40*10**3) * (Altitude_GB(iLon,iLat,iAlt,iBlock) - 130*10**3)
           endif 
           !See Gu et al., 2020 paper &
           !Fox and Hac, 2009 (Figure 6) for some really
           !good info on the branching ratios.
           avgHeatMatrix = MATMUL(MATMUL(Vibrational_Population, &
                Branching_Per_Vibration), Exotherm_O2_DR)
           avgHeat = avgHeatMatrix(1,1)
           
           ChemicalHeatingSub = &
               ChemicalHeatingSub + &
               Reaction * avgHeat

           !nu = 3-10 is a large fractional population ~0.2. How to handle heating?


           ! ----------------------------------------------------------
           ! O2 + hv ==> O + O 
           ! ----------------------------------------------------------

           Reaction = EuvDissRateS(iLon,iLat,iAlt,iO2_,iBlock)

           NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction
           NeutralSources(iO_) = NeutralSources(iO_) + 2.*Reaction

           !Dissociation heating is included in calc_euv.f90

           ! -----------------------------------------------------------
           ! CO2+ + e- ==> O + CO
              ! CO2+ + e- ==> CO(X^1 Sigma^+) + O(3P) + 8.31 eV  (a)
              ! CO2+ + e- ==> CO(a^3 Pi) + O(3P) + 2.30 eV        (b)
              ! CO2+ + e- ==> CO(a^'3 Sigma^+) + O(3P) + 1.26 eV (c)
              ! CO2+ + e- ==> CO(d^3 \Delta) + O(3P) + 0.49 eV     (d)
           !
           !According to Gu et al., 2020, originally (Rosati et al., 2003)
           !the appropriate branching ratios are
           !(a,b,c,d) = (0.24, 0.38, 0.18, 0.20)
           ! -----------------------------------------------------------
           
           
           Reaction =rtCO2P_e * Ions(ie_) * Ions(iCO2P_) 
          
           IonLosses(iCO2P_) = IonLosses(iCO2P_) + Reaction
           NeutralSources(iO_) = NeutralSources(iO_) + Reaction
           NeutralSources(iCO_) = NeutralSources(iCO_) + Reaction

           ChemicalHeatingSub = ChemicalHeatingSub + &
                0.24*Reaction*(8.31) + &
                0.38*Reaction*(2.30) + &
                0.18*Reaction*(1.26) + &
                0.20*Reaction*(0.49)

!-----------Ter-Molecular Neutral-Neutral Chemistry-------------------------------+
	   ! -----------------------------------------------------------
           ! O + CO + CO2 ==> 2.*CO2  
           ! -----------------------------------------------------------
           Reaction = rtO_CO_CO2 *Neutrals(iO_)*Neutrals(iCO_)*Neutrals(iCO2_)

           NeutralLosses(iO_) = NeutralLosses(iO_) + Reaction
           NeutralLosses(iCO_) = NeutralLosses(iCO_) + Reaction
           NeutralSources(iCO2_) = NeutralSources(iCO2_) + Reaction
           !reactionrate(9) = Reaction

	   ! -----------------------------------------------------------
           ! O + O + CO2 ==> O2 + CO2  
           ! -----------------------------------------------------------
           Reaction = rtO_O_CO2*Neutrals(iO_)*Neutrals(iO_)*Neutrals(iCO2_)

           NeutralLosses(iO_) = NeutralLosses(iO_) + 2.*Reaction
           NeutralSources(iO2_) = NeutralSources(iO2_) + Reaction

           ! No chemical heating included for O + O + M in Earth GITM

           ! -----------------------------------------------------------------------
           !     AIRGLOW QUANTITIES (SLAB WIDTH = 0.5*SCALE HEIGHT; cm units):
           !  *  O2 IR DAYGLOW FROM 3-BODY RECOMBINATION (via vtgcm2 coding)
           !  *  3-BODY O-RECOMBINATION INTO VARIOUS STATES: 
           !     (a) DIRECTLY to O2(a)-state and
           !     (b) INDIRECTLY to O2(5pi)-state plus quenching to O2(a)-state
           !     (c) INDIRECTLY to O2(c)-state plus quenching to O2(a)-state
           ! ---------------------------------------------------------------------
           
           !Volume Emission Rate computation (BP: 11/13/2020) from VTGCM
           !  ** VTGCM3 **  Fall 2007  (Slanger et al., 2006; Gerard et al., 2007)
           !  *  3-BODY O-RECOMBINATION (a), (b), (c). NET YIELD  Leading to O2(a)-state
           !     (Direct, Indirect leading to O2(a) state)
           ! beta = 0.75
           !  *  Quenching of O2(a)-state by CO2 and N2 
           quenchingFactor = 2.0e-20 * (Neutrals(iCO2_) + Neutrals(iN2_))
           fO2IR           = 1.0/(1.0 + 3800.0 * quenchingFactor)
           
           !Photons/cm3-sec
           o2VER           = Neutrals(iO_)*Neutrals(iO_)*Neutrals(iCO2_)*rtO_O_CO2* &
                             beta * fO2IR
           Emission(iEO2UV_) = Emission(iEO2UV_) + Reaction*1.0E-06 
!----------NO Nightglow------------------------------------------------------------+
           ! -----------------------------------------------------------
           ! N4S + O  ==> NO* + hv (Nightglow reaction)
           ! New code added by S. W. Bougher (170427)
           ! -----------------------------------------------------------
           Reaction = rtN4S_O*Neutrals(iN4S_)*Neutrals(iO_)

           NeutralLosses(iN4S_) = NeutralLosses(iN4S_) + Reaction
           NeutralLosses(iO_) = NeutralLosses(iO_) + Reaction
           NeutralSources(iNO_) = NeutralSources(iNO_) + Reaction

!          Total NO (delta plus gamma band system) mission (190-260 nm)
!          Initially Should be ph/m3.sec units  (check)
!          Convert to ph/cm3.sec units  
           Emission(iENOUV_) = Emission(iENOUV_) + Reaction*1.0E-06

           

           ! -----------------------------------------------------------
           ! NO + hv ==> N(4S) + O
           ! -----------------------------------------------------------

           Reaction = EuvDissRateS(iLon,iLat,iAlt,iNO_,iBlock) !rho * I * cross-sect
           
           NeutralLosses(iNO_) = NeutralLosses(iNO_) + Reaction
           NeutralSources(iN4S_) = NeutralSources(iN4S_) + Reaction
           NeutralSources(iO_) = NeutralSources(iO_) + Reaction

           !No chemical heating included in Earth GITM for this dissociation

           ! -----------------------------------------------------------
           ! N2D + O  ==> N4S + O
              ! -----------
              ! N(2D) + O -> N(4S) + O(3P) + 2.38 eV
              ! N(2D) + O -> N(4S) + O(1D) + 0.42 eV
              ! -----------
           ! -----------------------------------------------------------
           Reaction = rtN2D_O*Neutrals(iN2D_)*Neutrals(iO_)

           NeutralLosses(iN2D_) = NeutralLosses(iN2D_) + Reaction
           NeutralSources(iN4S_) = NeutralSources(iN4S_) + Reaction

           ChemicalHeatingSub = ChemicalHeatingSub + &
                                0.9 * Reaction * 2.38 + &
                                0.1 * Reaction * 0.42

           ! -----------------------------------------------------------
           ! N2D + CO2  ==> NO + CO
           ! -----------------------------------------------------------
           Reaction = rtN2D_CO2*Neutrals(iN2D_)*Neutrals(iCO2_)

           NeutralLosses(iN2D_) = NeutralLosses(iN2D_) + Reaction
           NeutralLosses(iCO2_) = NeutralLosses(iCO2_) + Reaction
!          NeutralLosses(iCO2_) = 0.0
           NeutralSources(iNO_) = NeutralSources(iNO_) + Reaction
           NeutralSources(iCO_) = NeutralSources(iCO_) + Reaction           

           ! -----------------------------------------------------------
           ! O2+ + N4S ==> NO+ + O
           ! -----------------------------------------------------------
           Reaction = rtO2P_N4S * Neutrals(iN4S_) * Ions(iO2P_)
           NeutralLosses(iN4S_) = NeutralLosses(iN4S_) + Reaction
           IonLosses(iO2P_) = IonLosses(iO2P_) + Reaction
           NeutralSources(iO_) = NeutralSources(iO_) + Reaction
           IonSources(iNOP_) = IonSources(iNOP_) + Reaction

            ! -----------------------------------------------------------
            ! NO+ + e- ==> O + N2D 
            ! g = 0.78
              ! -----------
              ! NO+ + e -> O + N(4S) + 2.75 eV  (0.22)
              ! NO+ + e -> O + N(2D) + 0.38 eV  (0.78)
              ! -----------
            ! -----------------------------------------------------------

            Reaction = rtNOP_e * Ions(ie_) * Ions(iNOP_) 
            IonLosses(iNOP_) = IonLosses(iNOP_) + Reaction
            NeutralSources(iN2D_) = NeutralSources(iN2D_) + Reaction
            NeutralSources(iO_) = NeutralSources(iO_) + Reaction

            ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * (0.38 * 0.78) + Reaction * (2.75 * 0.22)
            
            ! -----------------------------------------------------------
            ! NO+ + e- ==> O + N4S 
            ! (1-g) = 0.22
            ! -----------------------------------------------------------

            Reaction = rtNOP_e2 * Ions(ie_) * Ions(iNOP_) 
            IonLosses(iNOP_) = IonLosses(iNOP_) + Reaction
            NeutralSources(iN4S_) = NeutralSources(iN4S_) + Reaction
            NeutralSources(iO_) = NeutralSources(iO_) + Reaction 

            !Chemical heating taken care of in reaction above

!\
!-----------End Photochemistry-----------------------------------------------------+
!\

!\
!-----------Inert Species (zeroed out sources and sinks)
!-----------Species not calculated yet (zeroed out sources and sinks)
!/
              ! -----------------------------------------------------------
              ! He, H:  little idea about Aij coefficients
              ! -----------------------------------------------------------
!------------CO2 infinite reservoir assumption used-------------
!              NeutralLosses(iCO2_) = 0.0
!              NeutralSources(iCO2_) = 0.0
!------------Diffusion and advection only (no chemistry presently) -------
!              NeutralLosses(iAr_) = 0.0
!              NeutralSources(iAr_) = 0.0
!              NeutralLosses(iHe_) = 0.0
!              NeutralSources(iHe_) = 0.0
!              NeutralLosses(iH_) = 0.0
!              NeutralSources(iH_) = 0.0
!\
!------------------------------------------------------+
!-----------End Total Chemistry-------------------------------------------+
!------------------------------------------------------+

  !write(*,*) "CO2 balance - Sources, Losses"
  !write(*,*), NeutralSources(iCO2_), NeutralLosses(iCO2_)

  !write(*,*) "CO balance - Sources, Losses"
  !write(*,*), NeutralSources(iCO_), NeutralLosses(iCO_)

  !write(*,*) "O balance - Sources, Losses"
  !write(*,*), NeutralSources(iO_), NeutralLosses(iO_)

  !write(*,*) "O2 balance - Sources, Losses"
  !write(*,*), NeutralSources(iO2_), NeutralLosses(iO2_)  

if (.NOT. useNeutralChemistry) then 

   NeutralSources = 0.0
   NeutralLosses = 0.0
endif

call end_timing("chemsources")    
   end subroutine calc_chemical_sources



end Module ModChemistry

