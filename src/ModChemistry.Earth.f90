!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
 Module ModChemistry

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
    real :: tr, tr3, te3, te12d, te227d, te3m07, te12m056, te3m05, te3m039, te3m085
    real :: tr3m044,  tr3m08, tr3m04, te22m05, ti, tn, tn1, tn06, rr_opn2
    logical :: useNeutralConstituent(nSpeciesTotal),useIonConstituent(nIons)
    integer :: nittot = 0
contains
!-----------------------------------
subroutine calc_reaction_rates(iLon,iLat,iAlt,iBlock)
  
  integer, intent(in) :: iLon,iLat,iAlt,iBlock

   tn = Temperature(iLon,iLat,iAlt,iBlock)*&
       TempUnit(iLon,iLat,iAlt)
   
  tr = (iTemperature(iLon,iLat,iAlt,iBlock) &
       + tn) / 2.0

  tr3 = tr/300.0
  te3 = eTemperature(iLon,iLat,iAlt,iBlock)/300.0
  te12d = eTemperature(iLon,iLat,iAlt,iBlock)/1200.0
  te227d = -22740.0/eTemperature(iLon,iLat,iAlt,iBlock)

  te3m07   = te3**(-0.7)
  te12m056 = te12d**(-0.56)
  te3m05   = te3**(-0.5)
  te3m039  = te3**(-0.39)
  te3m085  = te3**(-0.85)
  
  tr3m044 = tr3**(-0.44)
  tr3m08  = tr3**(-0.8)
  tr3m04  = tr3**(-0.4)
  
  te22m05  = eTemperature(iLon,iLat,iAlt,iBlock)**(.5) * &
       exp(te227d)
  
  ti = iTemperature(iLon,iLat,iAlt,iBlock)
  
  
  tn1 = exp(107.8/tn)
  tn06 = exp(67.5/tn)
  
  rr_opn2 = min(5.0e-19,4.5e-20*tr3**2)
  

end subroutine calc_reaction_rates

subroutine calc_chemical_sources(iLon,iLat,iAlt,iBlock,IonSources, &
       IonLosses,NeutralSources, &
       NeutralLosses,ChemicalHeatingSub,Emission)

  use ModSources

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
 
  if (.not.UseIonChemistry) return

  IonSources = 0.0
  NeutralSources = 0.0
  IonLosses  = 0.0
  NeutralLosses = 0.0
  
  ChemicalHeatingSub = 0.0
  Emission = 0.0

  o2ptotal = 0
  o2total = 0
  ! ----------------------------------------------------------
  ! O2 -> 2O
  ! ----------------------------------------------------------

  rr=EuvDissRateS(iLon,iLat,iAlt,iO2_,iBlock)
  
  Reaction = rr * &
       Neutrals(iO2_)
  
  NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction
  NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + 2*Reaction

  ! -----------------------------------------------------------
  ! O+O+M -> O2+M
  ! -----------------------------------------------------------

  rr=9.59e-34 * exp(480./tn)
  rr= rr*1.e-12  !cm6s-1-->m6s-1
  
  Reaction = rr * Neutrals(iO_3P_)**2 *&
       (Neutrals(iO2_)+ &
       Neutrals(iO_3P_)+ &
       Neutrals(iN2_))

  NeutralLosses(iO_3P_) = NeutralLosses(iO_3P_) + 2*Reaction
  NeutralSources(iO2_) = NeutralSources(iO2_) + Reaction

  ! -----------------------------------------------------------
  ! O+O2+M -> O3+M
  ! -----------------------------------------------------------
!!$              rr=6.0e-34 * exp(300./tn)
!!$              rr=rr*1.e-12  !cm6s-1-->m6s-1
!!$
!!$              Reaction = rr * Neutrals(iO_3P_) *&
!!$                    Neutrals(iO2_) * &
!!$                   (Neutrals(iO2_) + &
!!$                    Neutrals(iO_3P_)  + &!$                    Neutrals(iN2_))
!!$
!!$              NeutralLosses(iO_3P_) = NeutralLosses(iO_3P_) + Reaction
!!$              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction
  

  ! ----------------------------------------------------------
  ! N2 -> 2N
  ! ----------------------------------------------------------
  
  Reaction = EuvDissRateS(iLon,iLat,iAlt,iN2_,iBlock) * &
       Neutrals(iN2_)
  
  !              Reaction = 5.0e-12 * &
  !                   Neutrals(iN2_)
  
  NeutralLosses(iN2_) = NeutralLosses(iN2_) + Reaction
  NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + .25*Reaction
  NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + .60*Reaction

  ! Solar EUV

  ! ----------------------------------------------------------
  ! N2+
  ! ----------------------------------------------------------
  
  ! Solar EUV
  
  Reaction = EuvIonRateS(iLon,iLat,iAlt,iN2P_,iBlock) * &
       Neutrals(iN2_)
  
  IonSources(iN2P_)   = IonSources(iN2P_)   + Reaction
  NeutralLosses(iN2_) = NeutralLosses(iN2_) + Reaction
  
  ! Aurora
  
  Reaction = AuroralIonRateS(iLon,iLat,iAlt,iN2_, iBlock) + &
       IonPrecipIonRateS(iLon,iLat,iAlt,iN2_, iBlock)
  
  IonSources(iN2P_)   = IonSources(iN2P_) + Reaction
  NeutralLosses(iN2_) = NeutralLosses(iN2_) + Reaction
  
  ! O+(2D) + N2 -> N2+ + O + 1.33 eV
  
  Reaction = &
       8.0e-16 * &
       Ions(iO_2DP_) * &
       Neutrals(iN2_)
  
  IonSources(iN2P_)    = IonSources(iN2P_)   + Reaction
  NeutralSources(iO_3P_)  = NeutralSources(iO_3P_) + Reaction
  IonLosses(iO_2DP_)   = IonLosses(iO_2DP_)  + Reaction
  NeutralLosses(iN2_)  = NeutralLosses(iN2_) + Reaction
  
  ChemicalHeatingSub = &
       ChemicalHeatingSub + Reaction * 1.33
  
  ChemicalHeatingS(iop2d_n2) =  &
        ChemicalHeatingS(iop2d_n2) + &
        Reaction * 1.33

   ! O+(2P) + N2 -> N2+ + O + 3.02 eV

  Reaction = &
       4.8e-16 * &
       Ions(iO_2PP_) * &
       Neutrals(iN2_)
  
  IonSources(iN2P_)    = IonSources(iN2P_)   + Reaction
  NeutralSources(iO_3P_)  = NeutralSources(iO_3P_) + Reaction
  IonLosses(iO_2PP_)   = IonLosses(iO_2PP_)  + Reaction
  NeutralLosses(iN2_)  = NeutralLosses(iN2_) + Reaction
  
  ChemicalHeatingSub = &
       ChemicalHeatingSub + Reaction * 3.02
  
  ChemicalHeatingS(iop2p_n2) =  &
       ChemicalHeatingS(iop2p_n2) + &
       Reaction * 3.02
  
  ! N2+ + O2 -> O2+ + N2 + 3.53 eV
  
  rr = 5.0e-17 * tr3m08
  
  Reaction = &
       rr * &
       Ions(iN2P_) * &
       Neutrals(iO2_)
  
  IonSources(iO2P_)    = IonSources(iO2P_)   + Reaction
  NeutralSources(iN2_) = NeutralSources(iN2_) + Reaction
  IonLosses(iN2P_)     = IonLosses(iN2P_)  + Reaction
  NeutralLosses(iO2_)  = NeutralLosses(iO2_) + Reaction
  
  o2ptotal = o2ptotal + reaction
  ChemicalHeatingSub = &
       ChemicalHeatingSub + Reaction * 3.53
  
  ChemicalHeatingS(in2p_o2) =  &
       ChemicalHeatingS(in2p_o2) + &
       Reaction * 3.53

  ! N2+ + O -> NO+ + N(2D) + 0.70 eV
!!!!!          -> NO+ + N(4S) + 3.08 eV
  
  rr = 1.4e-16 * tr3m044
  
  Reaction = &
       rr * &
       Ions(iN2P_) * &
       Neutrals(iO_3P_)
  
  IonSources(iNOP_)      = IonSources(iNOP_)      + Reaction
  NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + Reaction
  IonLosses(iN2P_)       = IonLosses(iN2P_)       + Reaction
  NeutralLosses(iO_3P_)     = NeutralLosses(iO_3P_)     + Reaction
  
  ChemicalHeatingSub = &
       ChemicalHeatingSub + &
       Reaction * 0.70 
  
  ChemicalHeatingS(in2p_o) =  &
       ChemicalHeatingS(in2p_o) + &
       Reaction * 0.70
  
  ! N2+ + e -> 2 N(2D) + 1.04 eV
  
  rr = 1.8e-13 * te3m039
  
  Reaction = &
       rr * &
       Ions(iN2P_) * &
       Ions(ie_)

  
  NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + 2*Reaction
  IonLosses(iN2P_)       = IonLosses(iN2P_)  + Reaction
  
  ChemicalHeatingSub = &
       ChemicalHeatingSub + &
       Reaction * 1.04
  
  ChemicalHeatingS(in2p_e) =  &
       ChemicalHeatingS(in2p_e) + &
       Reaction * 1.04
  
  ! N2+ + O -> O+(4S) + N2 + 1.96 eV
  
  rr = 1.4e-16 * tr3m044
  
  Reaction = &
       rr * &
       Ions(iN2P_) * &
       Neutrals(iO_3P_)
  
  NeutralSources(iN2_) = NeutralSources(iN2_) + Reaction
  IonSources(iO_4SP_)  = IonSources(iO_4SP_)  + Reaction
  NeutralLosses(iO_3P_)   = NeutralLosses(iO_3P_) + Reaction
  IonLosses(iN2P_)     = IonLosses(iN2P_)  + Reaction
  
  ChemicalHeatingSub = &
       ChemicalHeatingSub + &
       Reaction * 1.96
  
  ChemicalHeatingS(in2p_o) =  &
       ChemicalHeatingS(in2p_o) + &
       Reaction * 1.96
  
  ! N2+ + NO -> NO+ + N2 + 6.33 eV
  
  rr = 3.3e-16
  
  Reaction = &
       rr * &
       Ions(iN2P_) * &
       Neutrals(iNO_)
  
  NeutralSources(iN2_) = NeutralSources(iN2_) + Reaction
  IonSources(iNOP_)    = IonSources(iNOP_)     + Reaction
  NeutralLosses(iNO_)  = NeutralLosses(iNO_)  + Reaction
  IonLosses(iN2P_)     = IonLosses(iN2P_)    + Reaction
  
  ChemicalHeatingSub = &
       ChemicalHeatingSub + &
       Reaction * 6.33
  
  
  ! ----------------------------------------------------------
  ! O2+
  ! ----------------------------------------------------------
  
  ! -----------
  ! Solar EUV
  ! -----------
  
  Reaction = EuvIonRateS(iLon,iLat,iAlt,iO2P_,iBlock) * &
       Neutrals(iO2_)
  
  
  IonSources(iO2P_)   = IonSources(iO2P_)   + Reaction
  NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction
  
  o2ptotal = o2ptotal + reaction
  
  ! -----------
  ! Aurora
  ! -----------

  Reaction = AuroralIonRateS(iLon,iLat,iAlt,iO2_, iBlock) + &
       IonPrecipIonRateS(iLon,iLat,iAlt,iO2_, iBlock)
  
  IonSources(iO2P_)   = IonSources(iO2P_)   + Reaction
  
  o2ptotal = o2ptotal + reaction
  
  NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction


  ! -----------
  ! O+(4S) + O2 -> O2+ + O + 1.55 eV
  ! -----------
  
  rr = 2.0e-17 * tr3m04
  
  Reaction = &
       rr * &
       Ions(iO_4SP_) * &
       Neutrals(iO2_)
  
  NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
  IonSources(iO2P_)   = IonSources(iO2P_)   + Reaction
  o2ptotal = o2ptotal + reaction
  
  NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction
  IonLosses(iO_4SP_)  = IonLosses(iO_4SP_) + Reaction
  
  ChemicalHeatingSub = &
       ChemicalHeatingSub + &
       Reaction * 1.55
  
  ChemicalHeatingS(iop_o2) =  &
       ChemicalHeatingS(iop_o2) + &
       Reaction * 1.55
  
  ! -----------
  ! O+(2D) + O2 -> O2+ + O + 4.865 eV
  ! -----------
  
              rr = 7.0e-16

              Reaction = &
                   rr * &
                   Ions(iO_2DP_) * &
                   Neutrals(iO2_)

              NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
              IonSources(iO2P_)   = IonSources(iO2P_)   + Reaction
              o2ptotal = o2ptotal + reaction
              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction/1.0
              IonLosses(iO_2DP_)  = IonLosses(iO_2DP_) + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 4.865
              
              ! -----------
              ! N+ + O2 -> O2+ + N(4S) + 2.486 eV
              ! -----------

             rr = 1.1e-16

             Reaction = &
                  rr * &
                  Ions(iNP_) * &
                  Neutrals(iO2_)

             NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
             IonSources(iO2P_)      = IonSources(iO2P_)      + Reaction
              o2ptotal = o2ptotal + reaction
             NeutralLosses(iO2_)    = NeutralLosses(iO2_)    + Reaction
             IonLosses(iNP_)        = IonLosses(iNP_)       + Reaction

             ChemicalHeatingSub = &
                  ChemicalHeatingSub + &
                  Reaction * 2.486
             
              ChemicalHeatingS(inp_o2) =  &
                  ChemicalHeatingS(inp_o2) + &
                  Reaction * 2.486

              ! -----------
              ! N+ + O2 -> O2+ + N(2D) + 0.1 eV
              ! -----------

              rr = 2.0e-16

              Reaction = &
                   rr * &
                   Ions(iNP_) * &
                   Neutrals(iO2_)

             NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + Reaction

             IonSources(iO2P_)      = IonSources(iO2P_)      + Reaction
              o2ptotal = o2ptotal + reaction
             NeutralLosses(iO2_)    = NeutralLosses(iO2_)    + Reaction
             IonLosses(iNP_)        = IonLosses(iNP_)       + Reaction

             ChemicalHeatingSub = &
                  ChemicalHeatingSub + &
                  Reaction * 0.1

               ChemicalHeatingS(inp_o2) =  &
                  ChemicalHeatingS(inp_o2) + &
                  Reaction * 0.1


                 ! -----------
              ! O2+ + e -> O(1D) + O(1D) + 3.06 eV
              ! O2+ + e -> O(3P) + O(1D) + 5.02 eV
              ! O2+ + e -> O(3P) + O(3P) + 6.99 eV
              ! -----------

!              rr = 1.9e-13 * te3m05
              rr = 2.4e-13 * te3m07
!              if (eTemperature(iLon,iLat,iAlt,iBlock) .lt. 1200) then
!                 rr = 1.9e-13* te12m056
!              else
!                 rr = 7.4e-14* te12m056
!              endif

              Reaction = &
                   rr * &
                   Ions(iO2P_) * &
                   Ions(ie_)
!                   edensity
           
              if (UseNeutralConstituent(iO_1D_)) then
                 NeutralSources(iO_3P_)    = NeutralSources(iO_3P_) + 0.22*Reaction * 2.0
                 NeutralSources(iO_3P_)    = NeutralSources(iO_3P_) + 0.42*Reaction
                 NeutralSources(iO_1D_) = NeutralSources(iO_1D_) + 0.42*Reaction
                 NeutralSources(iO_1D_) = NeutralSources(iO_1D_) + 0.31*Reaction * 2.0

                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 6.99 * 0.22 + &
                      Reaction * 5.02 * 0.42 + &
                      Reaction * 3.06 * 0.31
                 
                 ChemicalHeatingS(io2p_e) = &
                      ChemicalHeatingS(io2p_e) + &
                      Reaction * 6.99 * 0.22 + &
                      Reaction * 5.02 * 0.42 + &
                      Reaction * 3.06 * 0.31
                 
              else
                 
                 NeutralSources(iO_3P_)    = NeutralSources(iO_3P_) + Reaction * 2.0

                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 5.0
                 
                  ChemicalHeatingS(io2p_e) = &
                      ChemicalHeatingS(io2p_e) + &
                      Reaction * 5.0

              endif
              
              IonLosses(iO2P_)      = IonLosses(iO2P_)   + Reaction

              ! -----------
              ! O2+ + N(4S) -> NO+ + O + 4.25 eV
              ! -----------

              rr = 1.8e-16

              Reaction = &
                   rr * &
                   Ions(iO2P_) * &
                   Neutrals(iN_4S_)

              NeutralSources(iO_3P_)   = NeutralSources(iO_3P_)   + Reaction
              IonSources(iNOP_)     = IonSources(iNOP_)     + Reaction
              NeutralLosses(iN_4S_) = NeutralLosses(iN_4S_) + Reaction
              IonLosses(iO2P_)      = IonLosses(iO2P_)      + Reaction
              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 4.25

              ChemicalHeatingS(io2p_n) =  &
                   ChemicalHeatingS(io2p_n) + &
                   Reaction * 4.25

              ! -----------
              ! O2+ + NO -> NO+ + O2 + 2.813 eV
              ! -----------

              rr = 4.4e-16

              Reaction = &
                   rr * &
                   Ions(iO2P_) * &
                   Neutrals(iNO_)

              NeutralSources(iO2_) = NeutralSources(iO2_) + Reaction
              IonSources(iNOP_)    = IonSources(iNOP_)    + Reaction
              NeutralLosses(iNO_)  = NeutralLosses(iNO_)  + Reaction
              IonLosses(iO2P_)     = IonLosses(iO2P_)     + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 2.813

              ChemicalHeatingS(io2p_no) =  &
                   ChemicalHeatingS(io2p_no) + &
                    Reaction * 2.813

              ! -----------
              ! O2+ + N2 -> NO+ + NO + 0.9333 eV
              ! -----------
              
              rr = 5.0e-22
              
              Reaction = &
                   rr * &
                   Ions(iO2P_) * &
                   Neutrals(iN2_)

              NeutralSources(iNO_) = NeutralSources(iNO_) + Reaction


              IonSources(iNOP_)    = IonSources(iNOP_)    + Reaction
              NeutralLosses(iN2_)  = NeutralLosses(iN2_)  + Reaction
              IonLosses(iO2P_)     = IonLosses(iO2P_)     + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 0.9333

               ChemicalHeatingS(io2p_n2) =  &
                   ChemicalHeatingS(io2p_n2) + &
                   Reaction * 0.9333


              ! ----------------------------------------------------------
              ! O(4S)+
              ! ----------------------------------------------------------

              ! Solar EUV

              Reaction = EuvIonRateS(iLon,iLat,iAlt,iO_4SP_,iBlock) * &
                   Neutrals(iO_3P_)

              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction

              ! Aurora

              Reaction = AuroralIonRateS(iLon,iLat,iAlt,iO_3P_, iBlock) + &
                  IonPrecipIonRateS(iLon,iLat,iAlt,iO_3P_, iBlock)

              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction

              ! -----------
              ! O+(2D) + O -> O+(4S) + O(3P) + 3.31 eV
              ! O+(2D) + O -> O+(4S) + O(1D) + 1.35 eV
              ! -----------

              rr = 1.0e-17

              Reaction = &
                   rr * &
                   Ions(iO_2DP_) * &
                   Neutrals(iO_3P_)

              ! We create and loose the same amount of O (when O(1D) is not used...
             
              if (UseNeutralConstituent(iO_1D_)) then
                 
                 NeutralSources(iO_3P_)  = NeutralSources(iO_3P_) + 0.5 * Reaction
                 NeutralSources(iO_1D_)  = NeutralSources(iO_1D_) + 0.5 * Reaction
                 
                 NeutralLosses(iO_3P_)   = NeutralSources(iO_3P_)  + Reaction

                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 3.31 * 0.5 + &
                      Reaction * 1.35 * 0.5
              else

                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 3.31

              endif
              
              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2DP_)  = IonLosses(iO_2DP_)  + Reaction
              
              ! -----------
              ! O+(2D) + e -> O+(4S) + e + 3.31 eV
              ! -----------

              rr = 7.8e-14 * te3m05

              Reaction = &
                   rr * &
                   Ions(iO_2DP_) * &
                   Ions(ie_)
!                   edensity

              ! We create and loose the same amount of e
              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2DP_)  = IonLosses(iO_2DP_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 3.31

              ChemicalHeatingS(iop2d_e) =  &
                   ChemicalHeatingS(iop2d_e) + &
                   Reaction * 3.31

              ! -----------
              ! O+(2D) + N2 -> O+(4S) + N2 + 3.31 eV
              ! -----------

              rr = 8.0e-16

              Reaction = &
                   rr * &
                   Ions(iO_2DP_) * &
                   Neutrals(iN2_)

              ! We create and loose the same amount of N2
              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2DP_)  = IonLosses(iO_2DP_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 3.31

              ChemicalHeatingS(iop2d_n2) =  &
                   ChemicalHeatingS(iop2d_n2) + &
                   Reaction * 3.31


              ! -----------
              ! O+(2P) + O -> O+(4S) + O + 5.0 eV
              ! -----------

              rr = 5.2e-17

              Reaction = &
                   rr * &
                   Ions(iO_2PP_) * &
                   Neutrals(iO_3P_)

              ! We create and loose the same amount of O
              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2PP_)  = IonLosses(iO_2PP_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 5.0

               ChemicalHeatingS(iop2p_o) =  &
                   ChemicalHeatingS(iop2p_o) + &
                    Reaction * 5.0
              ! -----------
              ! O+(2P) + e -> O+(4S) + e + 5.0 eV
              ! -----------

              rr = 4.0e-14 * te3m05

              Reaction = &
                   rr * &
                   Ions(iO_2PP_) * &
                   Ions(ie_)
!                   edensity

              ! We create and loose the same amount of e
              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2PP_)  = IonLosses(iO_2PP_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 5.0

              ChemicalHeatingS(iop2p_e) =  &
                   ChemicalHeatingS(iop2p_e) + &
                   Reaction * 5.0

              ! -----------
              ! O+(2P) -> O+(4S) + 2470A
              ! -----------

              rr = 0.047

              Reaction = &
                   rr * &
                   Ions(iO_2PP_)

              ! We create and loose the same amount of e
              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2PP_)  = IonLosses(iO_2PP_)  + Reaction

              Emission(iE2470_) = Emission(iE2470_) + Reaction

!              ! -----------
!              ! H+ + O -> O+(4S) + H + 0.0 eV
!              ! -----------
!
!              rr = 6.0e-10/1.0e6 * (8.0/9.0) * sqrt( &
!                   (ti + tn/4) * (tn + ti/16))
!
!              Reaction = &
!                   rr * &
!                   Ions(iHP_) * &
!                   Neutrals(iO_3P_)
!
!              NeutralSources(iH_) = NeutralSources(iH_) + Reaction
!              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
!              NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction
!              IonLosses(iHP_)     = IonLosses(iHP_)     + Reaction
!
!              ChemicalHeatingSub = &
!                   ChemicalHeatingSub + &
!                   Reaction * 0.0

              ! -----------
              ! N+ + O2 -> O+(4S) + NO + 2.31 eV
              ! -----------

              rr = 3.0e-17

              Reaction = &
                   rr * &
                   Ions(iNP_) * &
                   Neutrals(iO2_)

              NeutralSources(iNO_) = NeutralSources(iNO_) + Reaction
              IonSources(iO_4SP_)  = IonSources(iO_4SP_)  + Reaction
              NeutralLosses(iO2_)  = NeutralLosses(iO2_)  + Reaction
              IonLosses(iNP_)      = IonLosses(iNP_)      + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 2.31

               ChemicalHeatingS(inp_o2) =  &
                   ChemicalHeatingS(inp_o2) + &
                   Reaction * 2.31

              ! -----------
              ! O+(4S) + N2 -> NO+ + N(4S) + 1.10 eV
              ! -----------

              rr = rr_opn2

              Reaction = &
                   rr * &
                   Ions(iO_4SP_) * &
                   Neutrals(iN2_)

              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              IonSources(iNOP_)      = IonSources(iNOP_)      + Reaction
              NeutralLosses(iN2_)    = NeutralLosses(iN2_)    + Reaction
              IonLosses(iO_4SP_)     = IonLosses(iO_4SP_)     + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 1.10

              ChemicalHeatingS(iop_n2) =  &
                   ChemicalHeatingS(iop_n2) + &
                   Reaction * 1.10

              ! -----------
              ! O+(4S) + NO -> NO+ + O + 4.36 eV
              ! -----------

              rr = 8.0e-19

              Reaction = &
                   rr * &
                   Ions(iO_4SP_) * &
                   Neutrals(iNO_)

              NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
              IonSources(iNOP_)   = IonSources(iNOP_)   + Reaction
              NeutralLosses(iNO_) = NeutralLosses(iNO_) + Reaction
              IonLosses(iO_4SP_)  = IonLosses(iO_4SP_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 4.36

!              ! -----------
!              ! O+(4S) + H -> H+ + O + 0.0 eV
!              ! -----------
!
!              rr = 6.0e-10/1.0e6
!
!              Reaction = &
!                   rr * &
!                   Ions(iO_4SP_) * &
!                   Neutrals(iH_)
!
!              NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
!              IonSources(iHP_)    = IonSources(iHP_)    + Reaction
!              NeutralLosses(iH_)  = NeutralLosses(iH_)  + Reaction
!              IonLosses(iO_4SP_)  = IonLosses(iO_4SP_)  + Reaction
!
!              ChemicalHeatingSub = &
!                   ChemicalHeatingSub + &
!                   Reaction * 0.0

              ! -----------
              ! O+(4S) + N(2D) -> N+ + O + 1.45 eV
              ! -----------

              rr = 1.3e-16

              Reaction = &
                   rr * &
                   Ions(iO_4SP_) * &
                   Neutrals(iN_2D_)

              NeutralSources(iO_3P_)   = NeutralSources(iO_3P_)   + Reaction
              IonSources(iNP_)      = IonSources(iNP_)      + Reaction
              NeutralLosses(iN_2D_) = NeutralLosses(iN_2D_) + Reaction
              IonLosses(iO_4SP_)    = IonLosses(iO_4SP_)    + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 1.45
              
              ! ----------------------------------------------------------
              ! O(2D)+
              ! ----------------------------------------------------------

              ! Solar EUV

              Reaction = EuvIonRateS(iLon,iLat,iAlt,iO_2DP_,iBlock) * &
                   Neutrals(iO_3P_)

              IonSources(iO_2DP_) = IonSources(iO_2DP_) + Reaction
              NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction

              ! Aurora

              !Reaction = AuroralIonRateS(iLon,iLat,iAlt,iO_, iBlock)

              !IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              !NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction

              ! -----------
              ! O+(2P) + e -> O+(2D) + e + 1.69 eV
              ! -----------

              rr = 1.5e-13 * te3m05

              Reaction = &
                   rr * &
                   Ions(iO_2PP_) * &
                   Ions(ie_)
!                   edensity

              ! We create and loose the same amount of e
              IonSources(iO_2DP_) = IonSources(iO_2DP_) + Reaction
              IonLosses(iO_2PP_)  = IonLosses(iO_2PP_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 1.69

              ChemicalHeatingS(iop2p_e) =  &
                   ChemicalHeatingS(iop2p_e) + &
                   Reaction * 1.69
              
              ! -----------
              ! O+(2P) -> O+(2D) + 7320A
              ! -----------

              rr = 0.171

              Reaction = &
                   rr * &
                   Ions(iO_2PP_)

              IonSources(iO_2DP_) = IonSources(iO_2DP_) + Reaction
              IonLosses(iO_2PP_)  = IonLosses(iO_2PP_)  + Reaction

              Emission(iE7320_) = Emission(iE7320_) + Reaction

              ! -----------
              ! O+(2D) -> O+(4S) + 3726A
              ! -----------

              rr = 7.7e-5

              Reaction = &
                   rr * &
                   Ions(iO_2DP_)

              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2DP_)  = IonLosses(iO_2DP_)  + Reaction

              Emission(iE3726_) = Emission(iE3726_) + Reaction

              ! ----------------------------------------------------------
              ! O(2P)+
              ! ----------------------------------------------------------

              ! Solar EUV

              Reaction = EuvIonRateS(iLon,iLat,iAlt,iO_2PP_,iBlock) * &
                   Neutrals(iO_3P_)

              IonSources(iO_2PP_) = IonSources(iO_2PP_) + Reaction
              NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction

              ! Aurora

              !Reaction = AuroralIonRateS(iLon,iLat,iAlt,iO_, iBlock)

              !IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              !NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction

              ! -----------
              ! O+(2P) + N2 -> N+ + NO + 0.70 eV
              ! -----------

              rr = 1.0e-16

              Reaction = &
                   rr * &
                   Ions(iO_2PP_) * &
                   Neutrals(iN2_)

              NeutralSources(iNO_) = NeutralSources(iNO_) + Reaction
              IonSources(iNP_)     = IonSources(iNP_)     + Reaction
              NeutralLosses(iN2_)  = NeutralLosses(iN2_)  + Reaction
              IonLosses(iO_2PP_)   = IonLosses(iO_2PP_)   + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 0.70

               ChemicalHeatingS(iop2p_n2) =  &
                    ChemicalHeatingS(iop2p_n2) + &
                    Reaction * 0.70


              ! ----------------------------------------------------------
              ! N+
              ! ----------------------------------------------------------

              ! -----------
              ! O2+ + N(2D) -> N+ + O2 + 0.0 eV
              ! -----------

             rr = 2.5e-16

             Reaction = &
                  rr * &
                  Ions(iO2P_) * &
                  Neutrals(iN_2D_)

             NeutralSources(iO2_)  = NeutralSources(iO2_)  + Reaction
             IonSources(iNP_)      = IonSources(iNP_)      + Reaction
             NeutralLosses(iN_2D_) = NeutralLosses(iN_2D_) + Reaction
             IonLosses(iO2P_)      = IonLosses(iO2P_)      + Reaction
             ChemicalHeatingSub = &
                  ChemicalHeatingSub + &
                  Reaction * 0.0

!              ! -----------
!              ! He+ + N2 -> N+ + N + He + 0.28 eV
!              ! -----------
!
!              rr = 1.2e-9/1.0e6
!
!              Reaction = &
!                   rr * &
!                   Ions(iHeP_) * &
!                   Neutrals(iN2_)
!
!              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
!              NeutralSources(iHe_)   = NeutralSources(iHe_)   + Reaction
!              IonSources(iNP_)       = IonSources(iNP_)       + Reaction
!              NeutralLosses(iN2_)    = NeutralLosses(iN2_)    + Reaction
!              IonLosses(iHeP_)       = IonLosses(iHeP_)       + Reaction
!
!              ChemicalHeatingSub = &
!                   ChemicalHeatingSub + &
!                   Reaction * 0.28

              ! -----------
              ! O+(2P) + N -> N+ + O + 2.7 eV
              ! -----------

              rr = 1.0e-16

              Reaction = &
                   rr * &
                   Ions(iO_2PP_) * &
                   Neutrals(iN_4S_)

              NeutralSources(iO_3P_)   = NeutralSources(iO_3P_)   + Reaction
              IonSources(iNP_)      = IonSources(iNP_)      + Reaction
              NeutralLosses(iN_4S_) = NeutralLosses(iN_4S_) + Reaction
              IonLosses(iO_2PP_)    = IonLosses(iO_2PP_)    + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 2.7

              ! -----------
              ! O+(2D) + N -> N+ + O + 1.0 eV
              ! -----------

              rr = 7.5e-17

              Reaction = &
                   rr * &
                   Ions(iO_2DP_) * &
                   Neutrals(iN_4S_)

              NeutralSources(iO_3P_)   = NeutralSources(iO_3P_)   + Reaction
              IonSources(iNP_)      = IonSources(iNP_)      + Reaction
              NeutralLosses(iN_4S_) = NeutralLosses(iN_4S_) + Reaction
              IonLosses(iO_2DP_)    = IonLosses(iO_2DP_)    + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 1.0

              ! -----------
              ! N+ + O2 -> NO+ + O(1D) + 6.67 eV
              ! -----------

              rr = 2.6e-16

              Reaction = &
                   rr * &
                   Ions(iNP_) * &
                   Neutrals(iO2_)

              if (UseNeutralConstituent(iO_1D_)) then
                 
                 NeutralSources(iO_1D_) = NeutralSources(iO_1D_) + Reaction
                                                
              else

                 NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
               
              endif

              IonSources(iNOP_)   = IonSources(iNOP_)   + Reaction
              IonLosses(iNP_)     = IonLosses(iNP_)     + Reaction
              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 6.67

               ChemicalHeatingS(inp_o2) =  &
                   ChemicalHeatingS(inp_o2) + &
                   Reaction * 6.67

              ! -----------
              ! N+ + O -> O+ + N + 0.93 eV
              ! -----------

              rr = 5.0e-19

              Reaction = &
                   rr * &
                   Ions(iNP_) * &
                   Neutrals(iO_3P_)

              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              IonSources(iO_4SP_)    = IonSources(iO_4SP_)    + Reaction
              NeutralLosses(iO_3P_)     = NeutralLosses(iO_3P_)     + Reaction
              IonLosses(iNP_)        = IonLosses(iNP_)        + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 0.93

              ChemicalHeatingS(inp_o) =  &
                   ChemicalHeatingS(inp_o) + &
                   Reaction * 0.93

!              ! -----------
!              ! N+ + H -> H+ + N + 0.90 eV
!              ! -----------
!
!              rr = 3.6e-12/1.0e6
!
!              Reaction = &
!                   rr * &
!                   Ions(iNP_) * &
!                   Neutrals(iH_)
!
!              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
!              IonSources(iHP_)       = IonSources(iHP_)       + Reaction
!              NeutralLosses(iH_)     = NeutralLosses(iH_)     + Reaction
!              IonLosses(iNP_)        = IonLosses(iNP_)        + Reaction
!
!              ChemicalHeatingSub = &
!                   ChemicalHeatingSub + &
!                   Reaction * 0.9

              ! ----------------------------------------------------------
              ! NO+
              ! ----------------------------------------------------------

              ! -----------
              !!!!!!! NO+ + e -> O + N(4S) + 2.75 eV  (0.22)
              ! NO+ + e -> O + N(2D) + 0.38 eV  (0.78)
              ! -----------

              rr = 4.2e-13 * te3m085

              Reaction = &
                   rr * &
                   Ions(iNOP_) * &
                   Ions(ie_)
!                   edensity


              NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + Reaction
              NeutralSources(iO_3P_)    = NeutralSources(iO_3P_)    +      Reaction
              IonLosses(iNOP_)       = IonLosses(iNOP_)       +      Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * (0.38 * 0.78)

               ChemicalHeatingS(inop_e) =  &
                   ChemicalHeatingS(inop_e) + &
                   Reaction * (0.38 * 0.78)

              ! ----------------------------------------------------------
              ! N(4S)
              ! ----------------------------------------------------------

              ! -----------
              ! N(2D) + e -> N(4S) + e + 2.38 eV
              ! -----------

              rr = 5.5e-16 * te3 ** (0.5)

              Reaction = &
                   rr * &
                   Neutrals(iN_2D_) * &
                   Ions(ie_)
!                   edensity

              ! We create and loose the same amount of e
              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              NeutralLosses(iN_2D_)  = NeutralLosses(iN_2D_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 2.38

              ChemicalHeatingS(in2d_e) =  &
                   ChemicalHeatingS(in2d_e) + &
                   Reaction * 2.38

              ! -----------
              ! N(2D) + O -> N(4S) + O(3P) + 2.38 eV
              ! N(2D) + O -> N(4S) + O(1D) + 0.42 eV
              ! -----------

              rr = 2.0e-18

              Reaction = &
                   rr * &
                   Neutrals(iN_2D_) * &
                   Neutrals(iO_3P_)

              if (UseNeutralConstituent(iO_1D_)) then
                 NeutralSources(iO_3P_)    = NeutralSources(iO_3P_)     + 0.9 * Reaction
                 NeutralSources(iO_1D_)    = NeutralSources(iO_1D_)  + 0.1 * Reaction
                 NeutralLosses(iO_3P_)     = NeutralLosses(iO_3P_)  + Reaction
                 
                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      0.9 * Reaction * 2.38 + &
                      0.1 * Reaction * 0.42
                 
                 ChemicalHeatingS(in2d_o) =  &
                      ChemicalHeatingS(in2d_o) + &
                      0.9 * Reaction * 2.38 + &
                      0.1 * Reaction * 0.42

              else
                 
                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 2.38

              endif

              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              NeutralLosses(iN_2D_)  = NeutralLosses(iN_2D_)  + Reaction
              
              ! -----------
              ! N(2D) -> N(4S) + 5200A
              ! -----------

              rr = 1.06e-5

              Reaction = &
                   rr * &
                   Neutrals(iN_2D_)

              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              NeutralLosses(iN_2D_)  = NeutralLosses(iN_2D_)  + Reaction
              Emission(iE5200_) = Emission(iE5200_) + Reaction

              ! -----------
              ! NO -> N(4S) + O
              ! -----------

!              rr = 8.3e-6
              rr=4.5e-6*exp(-1.e-8*(Neutrals(iO2_)*1.e-6)**0.38)

              Reaction = &
                   rr * &
                   Neutrals(iNO_)

              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              NeutralSources(iO_3P_)    = NeutralSources(iO_3P_)    + Reaction
              NeutralLosses(iNO_)    = NeutralLosses(iNO_)    + Reaction

              ! -----------
              ! N(4S) + O2 -> NO + O + 1.385 eV
              ! -----------
              
              rr = 4.4e-18 * exp(-3220/tn) 

              Reaction = &
                   rr * &
                   Neutrals(iN_4S_) * &
                   Neutrals(iO2_)

              NeutralSources(iNO_)  = NeutralSources(iNO_)  + Reaction
              NeutralSources(iO_3P_)   = NeutralSources(iO_3P_)   + Reaction
              NeutralLosses(iN_4S_) = NeutralLosses(iN_4S_) + Reaction
              NeutralLosses(iO2_)   = NeutralLosses(iO2_)   + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 1.385

              ChemicalHeatingS(in_o2) =  &
                   ChemicalHeatingS(in_o2) + &
                   Reaction * 1.385

              ! -----------
              ! N(4S) + NO -> N2 + O + 3.25 eV
              ! -----------
              rr = 1.5e-18 * sqrt(tn)
!              rr = 3.4e-17


              Reaction = &
                   rr * &
                   Neutrals(iN_4S_) * &
                   Neutrals(iNO_)

              NeutralSources(iN2_)  = NeutralSources(iN2_)  + Reaction
              NeutralSources(iO_3P_)   = NeutralSources(iO_3P_)   + Reaction
              NeutralLosses(iN_4S_) = NeutralLosses(iN_4S_) + Reaction
              NeutralLosses(iNO_)   = NeutralLosses(iNO_)   + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 3.25
              
               ChemicalHeatingS(ino_n) =  &
                   ChemicalHeatingS(ino_n) + &
                    Reaction * 3.25


              ! -----------
              ! N(2P) -> N(2D) + 10400A
              ! -----------

              rr = 7.9e-2

              Reaction = &
                   rr * &
                   Neutrals(iN_2P_)

              NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + Reaction
              NeutralLosses(iN_2P_)  = NeutralLosses(iN_2P_)  + Reaction

              Emission(iE10400_) = Emission(iE10400_) + Reaction

              ! ----------------------------------------------------------
              ! N(2D)
              ! ----------------------------------------------------------

              ! -----------
              ! N(2D) + O2 -> NO + O(3P) + 3.76 eV
              ! N(2D) + O2 -> NO + O(1D) + 1.80 eV
              ! -----------

              rr = 6.2e-18 *(Tn/300)

              Reaction = &
                   rr * &
                   Neutrals(iN_2D_) * &
                   Neutrals(iO2_)

              if (UseNeutralConstituent(iO_1D_)) then
                 
                 NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + 0.9 * Reaction
                 NeutralSources(iO_1D_) = NeutralSources(iO_1D_) + 0.1 * Reaction
                 
                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 3.76 * 0.9 + &
                      Reaction * 1.80 * 0.1
                 
                 ChemicalHeatingS(in2d_o2) =  &
                      ChemicalHeatingS(in2d_o2) + &
                      Reaction * 3.76 * 0.9 + &
                      Reaction * 1.80 * 0.1

              else
                 
                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 3.76
                 
              endif
              NeutralSources(iNO_)   = NeutralSources(iNO_)   + Reaction
              NeutralLosses(iN_2D_) = NeutralLosses(iN_2D_) + Reaction
              NeutralLosses(iO2_)   = NeutralLosses(iO2_)   + Reaction

              ! -----------
              ! N(2D) + NO -> N2 + O + 5.63 eV
              ! -----------

              rr = 7.0e-17

              Reaction = &
                   rr * &
                   Neutrals(iN_2D_) * &
                   Neutrals(iNO_)

              NeutralSources(iN2_)  = NeutralSources(iN2_)  + Reaction
              NeutralSources(iO_3P_)   = NeutralSources(iO_3P_)   + Reaction
              NeutralLosses(iN_2D_) = NeutralLosses(iN_2D_) + Reaction
              NeutralLosses(iNO_)   = NeutralLosses(iNO_)   + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 5.63

              ! ----------------------------------------------------------
              ! O(1D)
              ! ----------------------------------------------------------

              if (UseNeutralConstituent(iO_1D_)) then
                 ! ------------
                 ! O(1D) -> O(3P) + 6300A
                 ! ------------
                 
                 rr = 0.0071
                 
                 Reaction = &
                      rr * &
                      Neutrals(iO_1D_)
                 
                 NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
                 NeutralLosses(iO_1D_)  = NeutralLosses(iO_1D_)  + Reaction
                 
                 Emission(iE6300_) = Emission(iE6300_) + Reaction              
                 
                 
                 ! ------------                                                             
                 ! O(1D) -> O(3P) + 6364A                                               
                 ! ------------                                                             
                 
                 rr = 0.0022
                 
                 Reaction = &
                      rr * &
                      Neutrals(iO_1D_)
                 
                 NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
                 NeutralLosses(iO_1D_)  = NeutralLosses(iO_1D_)  + Reaction
                 
                 Emission(iE6364_) = Emission(iE6364_) + Reaction
                 
                 ! ------------                                                            
                 ! O(1D) + e -> O(3P) + e
                 ! ------------                                                            
                 
                 rr = 2.6e-17 * te22m05
                 
                 Reaction = &
                      rr * &
                      Neutrals(iO_1D_)
                 
                 NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
                 NeutralLosses(iO_1D_)  = NeutralLosses(iO_1D_)  + Reaction
                 
                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 1.96
                 
                 ChemicalHeatingS(io1d_e) = &
                      ChemicalHeatingS(io1d_e) + &
                      Reaction * 1.96
                 
                 ! ------------
                 ! O(1D) + N2 -> O(3P) + N2 + 1.96 eV
                 ! ------------
                 
                 rr = 2.3e-17
                 Reaction = &
                      rr * &
                      Neutrals(iO_1D_) * &
                      Neutrals(iN2_)
                 
                 !              !We create and loose the same amount of N2
                 NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
                 NeutralLosses(iO_1D_) = NeutralLosses(iO_1D_) + Reaction
                 
                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 1.96
                 
                 ChemicalHeatingS(io1d_n2) = &
                      ChemicalHeatingS(io1d_n2) + &
                      Reaction * 1.96
                 
                 
                 ! ------------
                 ! O(1D) + O2 -> O(3P) + O2 + 1.96 eV
                 ! ------------
                 
                 rr = 2.9e-17 * tn06
                 
                 Reaction = &
                      rr * &
                      Neutrals(iO_1D_) * &
                      Neutrals(iO2_)
                 
                 !We create and loose the same amount of O2                              
                 NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
                 NeutralLosses(iO_1D_) = NeutralLosses(iO_1D_) + Reaction
                 
                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 1.96
                 
                 ChemicalHeatingS(io1d_o2) = &
                      ChemicalHeatingS(io1d_o2) + &
                      Reaction * 1.96
                 
                 ! ------------                                                           
                 ! O(1D) + O(3P) -> O(3P) + O(3P) + 1.96 eV
                 ! ------------                                                            
                 
                 rr = 8.0e-18
                 Reaction = &
                      rr * &
                      Neutrals(iO_1D_) * &
                      Neutrals(iO_3P_)
                 
                 NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
                 NeutralLosses(iO_1D_) = NeutralLosses(iO_1D_) + Reaction
                 
                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 1.96
                 
                 ChemicalHeatingS(io1d_o) = &
                      ChemicalHeatingS(io1d_o) + &
                      Reaction * 1.96
                 
              endif
              
              ! ----------------------------------------------------------
              ! NO
              ! ----------------------------------------------------------
              ! -----------
              ! NO -> NO+ + e
              ! -----------

!              rr = 6.0e-7

              rr=5.88e-7*(1+0.2*(f107-65)/100)*exp(-2.115e-18* &
                   (Neutrals(iO2_)*1.e-6)**0.8855)

              Reaction = &
                   rr * &
                   Neutrals(iNO_)

              IonSources(iNOP_)   = IonSources(iNOP_)   + Reaction
              NeutralLosses(iNO_) = NeutralLosses(iNO_) + Reaction


    
   end subroutine calc_chemical_sources



end Module ModChemistry

