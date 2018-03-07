!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine calc_chemistry(iBlock)

  use ModSizeGitm
  use ModGITM
  use ModPlanet
  use ModRates
  use ModEUV
  use ModSources
  use ModInputs, only: iDebugLevel, UseIonChemistry, UseNeutralChemistry,f107,f107a
  use ModConstants
  use ModTime, only: istep,utime,ijulianday
  use EUA_ModMsis90, only: meter6, gtd6
  implicit none

  integer, intent(in) :: iBlock

  real :: IonSources(nIons), NeutralSources(nSpeciesTotal)
  real :: IonLosses(nIons), NeutralLosses(nSpeciesTotal)
  real :: DtSub, DtOld, DtTotal, DtMin, DtAve, Source, Reaction, tr, tr3, rr
  real :: te3, ti, tn, tn1, tn06, dtsubtmp, losstmp, dentmp, l, t,m1,m2,y1,y2,k1,k2
  real :: Ions(nIons), Neutrals(nSpeciesTotal)
  real :: tli(nIons), tsi(nIons), tln(nSpeciesTotal), tsn(nSpeciesTotal)
  real :: szap

  integer :: iLon, iLat, iAlt, iIon, nIters, iNeutral
  
  real :: lon, lat, alt
  real, dimension(1:2) :: msis_temp
  real, dimension(1:8) :: msis_dens
  real :: LonDeg, LatDeg, AltKm, LST
  real,dimension(7)    :: AP  

  real :: ChemicalHeatingSub,percent,o2ptotal
  real :: ChemicalHeatingSubI, ChemicalHeatingSubE
  real :: Emission(nEmissions), EmissionTotal(nEmissions)
  real, dimension(nLons,nLats,nAlts) :: &
       tr3d, tr33d,te12d, tr3m0443d, tr3m083d, tr3m043d, &
       tn3d, ti3d, ti33d, ti103d,ti93d, ti153d, ti3m0443d, ti3m0243d, &
       ti3m0233d, ti3m1163d, ti10m0673d, ti3m0453d, ti10m2123d, ti3m0873d, &
       ti3m0523d, ti9m0923d, ti3m0553d, ti15m0233d, &
       te3m0393d, te3m0853d, te33d, te3m053d,te3m073d,&
       te12m0563d,te227d, te3m0913d, te3m0813d, te073d

  real, dimension(nLons,nLats,nAlts) :: &
       teffective_n2, teffective_o2, teffective_no, u2, mb, mbb, &
       k1_n2, k2_o2, k3_no

  real :: k1_n2_point, k2_o2_point, k3_no_point

  real :: te3m05,te3m07,te12m056, tr3m044, tr3m04, tr3m08, te07
  real :: ti3m044, ti3m024, ti3m023, ti3m116, ti10m067, ti3m045
  real :: ti10m212, ti3m087, ti3m052, ti9m092, ti3m055, ti15m023
  real :: te3m039, te3m085, rr_opn2, te22m05, te3m091, te3m081
  real :: ionso, ionlo, neuso, neulo

  logical :: UseNeutralConstituent(nSpeciesTotal)
  logical :: UseIonConstituent(nIons)
  !---------------------------------------------------------------------------

  if (iDebugLevel > 3) then
     do iIon = 1, nIons
        write(*,*) "====> start calc_chemistry: Max Ion Density: ", iIon, &
             maxval(IDensityS(1:nLons,1:nLats,(nAlts*4)/5,iIon,iBlock))
     enddo
  endif

 
  UseNeutralConstituent = .true.
  UseIonConstituent     = .true.

!  UseNeutralConstituent(iO_1D_) = .false.
!  UseIonConstituent(iO_2PP_) = .false.
!  UseIonConstituent(iO_2DP_) = .false.
!  
!  UseNeutralConstituent(iN_4S_) = .false.
!  UseNeutralConstituent(iN_2D_) = .false.
!  UseNeutralConstituent(iO2_) = .false.
!
!  UseNeutralConstituent = .false.
!  UseIonConstituent = .false.
!  UseIonConstituent(1) = .true.
!  UseIonConstituent(2) = .true.
!  UseIonConstituent(3) = .true.
!
!  UseNeutralConstituent(1) = .true.
!  UseNeutralConstituent(2) = .true.
!  UseNeutralConstituent(3) = .true.

!  open(unit=95,file='data.dat')
  DtMin = Dt

  if (.not.UseIonChemistry) return

  call report("Chemistry",2)
  call start_timing("calc_chemistry")

  DtAve = 0.0

  nIters=0

!  AuroralIonRateS = 0.0

  u2 = IVelocity(1:nLons,1:nLats,1:nAlts,iEast_,iBlock)**2 + &
       IVelocity(1:nLons,1:nLats,1:nAlts,iNorth_,iBlock)**2 + &
       IVelocity(1:nLons,1:nLats,1:nAlts,iUp_,iBlock)**2

  mb  = 0.0
  mbb = 0.0

  ! This is from Schunk and Nagy 2nd Ed, formula 12.13 (pg 416)
  do iNeutral = 1, nSpeciesTotal
     ! Collisions should be better defined
     mbb = mbb + &
          (Collisions(1:nLons,1:nLats,1:nAlts,iVIN_)) / &
          (mass(iNeutral) + MassI(iO_4SP_))
     mb  = mb + &
          (mass(iNeutral) * Collisions(1:nLons,1:nLats,1:nAlts,iVIN_)) / &
          (mass(iNeutral) + MassI(iO_4SP_))
  enddo

  mb = mb/mbb

  teffective_n2 = iTemperature(1:nLons,1:nLats,1:nAlts,iBlock) + &
       MassI(iO_4SP_)/(MassI(iO_4SP_) - Mass(iN2_)) * &
       (Mass(iN2_) + mb)/(3*Boltzmanns_Constant) * u2
       
  teffective_o2 = iTemperature(1:nLons,1:nLats,1:nAlts,iBlock) + &
       MassI(iO_4SP_)/(MassI(iO_4SP_) - Mass(iO2_)) * &
       (Mass(iO2_) + mb)/(3*Boltzmanns_Constant) * u2
       
  teffective_no = iTemperature(1:nLons,1:nLats,1:nAlts,iBlock) + &
       MassI(iO_4SP_)/(MassI(iO_4SP_) - Mass(iNO_)) * &
       (Mass(iNO_) + mb)/(3*Boltzmanns_Constant) * u2
       
  where (teffective_n2 < 350.0)
     teffective_n2 = 350.0
  endwhere

  where (teffective_n2 > 6000.0)
     teffective_n2 = 6000.0
  endwhere

  where (teffective_o2 < 350.0)
     teffective_o2 = 350.0
  endwhere

  where (teffective_o2 > 6000.0)
     teffective_o2 = 6000.0
  endwhere

  where (teffective_no < 320.0)
     teffective_no = 320.0
  endwhere

  where (teffective_no > 6000.0)
     teffective_no = 6000.0
  endwhere

  where (teffective_n2 <= 1700.0)
     k1_n2 =   1.533e-18 &
          - 5.920e-19*(teffective_n2/300.0) &
          + 8.600e-20*(teffective_n2/300.0)**2
  endwhere
  
  where (teffective_n2 > 1700.0)
     k1_n2 =   2.730e-18 &
          - 1.155e-18*(teffective_n2/300.0) &
          + 1.483e-19*(teffective_n2/300.0)**2
  endwhere
  
  k2_o2 =   2.820e-17 &
       - 7.740e-18*(teffective_o2/300.0) &
       + 1.073e-18*(teffective_o2/300.0)**2 &
       - 5.170e-20*(teffective_o2/300.0)**3 &
       + 9.650e-22*(teffective_o2/300.0)**4
  
  where (teffective_no <= 1500.0)
     k3_no =   8.360e-19 &
          - 2.020e-19*(teffective_no/300.0) &
          + 6.950e-20*(teffective_no/300.0)**2
  endwhere
  
  where (teffective_no > 1500.0)
     k3_no =   5.330e-19 &
          - 1.640e-20*(teffective_no/300.0) &
          + 4.720e-20*(teffective_no/300.0)**2 &
          - 7.050e-22*(teffective_no/300.0)**3
  endwhere
  
   tr3d = (iTemperature(1:nLons,1:nLats,1:nAlts,iBlock) &
          + Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*&
          TempUnit(1:nLons,1:nLats,1:nAlts)) / 2.0

  tn3d = Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*&
                TempUnit(1:nLons,1:nLats,1:nAlts)
				
  ti3d = 	iTemperature(1:nLons,1:nLats,1:nAlts,iBlock)			
		 
		 
  ti33d = ti3d/300.0
  ti103d = ti3d/1000.0
  ti93d = ti3d/900.0
  ti153d = ti3d/1500.0
  te33d = eTemperature(1:nLons,1:nLats,1:nAlts,iBlock)/300.0
  te12d = eTemperature(1:nLons,1:nLats,1:nAlts,iBlock)/1200.0
  te227d = -22740.0/eTemperature(1:nLons,1:nLats,1:nAlts,iBlock)
  te073d = (250.0/eTemperature(1:nLons,1:nLats,1:nAlts,iBlock))**0.7

  te3m073d    = te33d**(-0.7)
  te12m0563d  = te12d**(-0.56)
  te3m053d    = te33d**(-0.5)
  te3m0393d   = te33d**(-0.39)
  te3m0853d   = te33d**(-0.85)
  te3m0913d   = te33d**(0.91)
  te3m0813d   = te33d**(0.81)
  ti3m0443d   = ti33d**(-0.44)
  ti3m0243d   = ti33d**(-0.24)
  ti3m0233d   = ti33d**(-0.23)
  ti3m1163d   = ti33d**(-1.16)
  ti10m0673d  = ti103d**(0.67)
  ti3m0453d   = ti33d**(-0.45)
  ti10m2123d  = ti103d**(2.12)
  ti3m0873d   = ti33d**(0.87)
  ti3m0523d   = ti33d**(-0.52)
  ti9m0923d   = ti93d**(0.92)
  ti3m0553d   = ti33d**(0.55)
  ti15m0233d  = ti153d**(0.2)

  m1 = ALOG(1.0/100000.0)/(115.0-100.0)
  k1 = 100000.0*exp(-m1*100.0)
  m2 = ALOG(1.0/1000.0)/(180.0-100.0)
  k2 = 1000.0*exp(-m2*100.0)


  do iLon = 1, nLons
     do iLat = 1, nLats

        szap = cos(sza(iLon, iLat,iBlock))
        if (szap < 0.0) szap = 0.0

        do iAlt = 1, nAlts

           y1 = max(1.0,k1*exp(m1*altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0))
           y2 = max(1.0,k2*exp(m2*altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0))
           NeutralSourcesTotal = 0.0
           NeutralLossesTotal = 0.0
           tr       = tr3d(iLon,iLat,iAlt)
           tr3      = tr33d(iLon,iLat,iAlt)
           te3      = te33d(iLon,iLat,iAlt)
           te3m05   = te3m053d(iLon,iLat,iAlt)
           te3m07   = te3m073d(iLon,iLat,iAlt)
           te07     = te073d(iLon,iLat,iAlt)
           te3m085  = te3m0853d(iLon,iLat,iAlt)
           te3m091  = te3m0913d(iLon,iLat,iAlt)
           te3m081  = te3m0813d(iLon,iLat,iAlt)
           te3m039  = te3m0393d(iLon,iLat,iAlt)
           te12m056 = te12m0563d(iLon,iLat,iAlt)
           tr3m044  = tr3m0443d(iLon,iLat,iAlt)
           tr3m04   = tr3m043d(iLon,iLat,iAlt)
           ti3m044  = ti3m0443d (iLon,iLat,iAlt)
           ti3m024  = ti3m0243d (iLon,iLat,iAlt)
           ti3m023  = ti3m0233d (iLon,iLat,iAlt)
           ti3m116  = ti3m1163d (iLon,iLat,iAlt)
           ti10m067 = ti10m0673d (iLon,iLat,iAlt)
           ti3m045  = ti3m0453d (iLon,iLat,iAlt)
           ti10m212 = ti10m2123d (iLon,iLat,iAlt)
           ti3m087  = ti3m0873d (iLon,iLat,iAlt)
           ti3m052  = ti3m0523d (iLon,iLat,iAlt)
           ti9m092  = ti9m0923d (iLon,iLat,iAlt)
           ti3m055  = ti3m0553d (iLon,iLat,iAlt)
           ti15m023 = ti15m0233d (iLon,iLat,iAlt)
           te22m05  = eTemperature(iLon,iLat,iAlt,iBlock)**(.5) * &
                exp(te227d(iLon,iLat,iAlt))
           tr3m08  = tr3m083d(iLon,iLat,iAlt)
           ti = iTemperature(iLon,iLat,iAlt,iBlock)
           tn = Temperature(iLon,iLat,iAlt,iBlock)*&
                TempUnit(iLon,iLat,iAlt)

           tn1 = exp(107.8/tn)
           tn06 = exp(67.5/tn)

           rr_opn2 = min(5.0e-19,4.5e-20*tr3**2)

           k1_n2_point = k1_n2(iLon,iLat,iAlt)
           k2_o2_point = k2_o2(iLon,iLat,iAlt)
           k3_no_point = k3_no(iLon,iLat,iAlt)

           DtTotal = 0.0
           EmissionTotal = 0.0

           Ions = IDensityS(iLon,iLat,iAlt,:,iBlock)

           Neutrals = NDensityS(iLon,iLat,iAlt,:,iBlock)

           niters = 0
           o2ptotal = 0

           do while (DtTotal < Dt)
              
              ChemicalHeatingSub = 0.0
              ChemicalHeatingSubI = 0.0
              ChemicalHeatingSubE = 0.0
              ChemicalHeatingS = 0
              Emission = 0.0

              DtSub = Dt - DtTotal

              IonSources = 0.0
              NeutralSources = 0.0
              IonLosses  = 0.0
              NeutralLosses = 0.0

              ! ----------------------------------------------------------
              ! O2 -> 2O
              ! ----------------------------------------------------------
              rr=EuvDissRateS(iLon,iLat,iAlt,iO2_,iBlock) + &
                  O2PERateS(iLon,iLat,iAlt,1,iBlock)

              Reaction = rr * &
                         Neutrals(iO2_)

              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction
              NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + 2*Reaction

              ! -----------------------------------------------------------
              ! O+O+M -> O2+M + 5.12 eV
              ! -----------------------------------------------------------
              rr=4.7e-33 *((300./tn)**(2))
              rr= rr*1.e-12  !cm6s-1-->m6s-1

              Reaction = rr * Neutrals(iO_3P_)**2 *&
                   (Neutrals(iO2_)+ &
                    Neutrals(iO_3P_)+ &
                    Neutrals(iN2_))

              NeutralLosses(iO_3P_) = NeutralLosses(iO_3P_) + 2*Reaction
              NeutralSources(iO2_) = NeutralSources(iO2_) + Reaction
			  
              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 5.12 

              ! ----------------------------------------------------------
              ! N2 -> 2N
              ! ----------------------------------------------------------

              Reaction = EuvDissRateS(iLon,iLat,iAlt,iN2_,iBlock) * &
                          Neutrals(iN2_) + &
                          N2PERateS(iLon,iLat,iAlt,1,iBlock)*Neutrals(iN2_)
              !rr=N2DissRateS(iLon,iLat,iAlt,1,iBlock) + &
              !     N2PERateS(iLon,iLat,iAlt,1,iBlock)

              NeutralLosses(iN2_) = NeutralLosses(iN2_) + Reaction
              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + .25*Reaction
              NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + .60*Reaction

              ! Solar EUV

              ! ----------------------------------------------------------
              ! N2+
              ! ----------------------------------------------------------

              Reaction = EuvIonRateS(iLon,iLat,iAlt,iN2P_,iBlock) * &
                          Neutrals(iN2_) + &
                          N2PERateS(iLon,iLat,iAlt,3,iBlock)*Neutrals(iN2_)

              IonSources(iN2P_)   = IonSources(iN2P_)   + Reaction
              NeutralLosses(iN2_) = NeutralLosses(iN2_) + Reaction

              ! Aurora

              Reaction = AuroralIonRateS(iLon,iLat,iAlt,iN2_, iBlock) + &
                   IonPrecipIonRateS(iLon,iLat,iAlt,iN2_, iBlock)

              IonSources(iN2P_)   = IonSources(iN2P_) + Reaction
              NeutralLosses(iN2_) = NeutralLosses(iN2_) + Reaction

              ! O+(2D) + N2 -> N2+ + O + 1.35 eV
			  
			  rr = 1.5e-16 * ti3m055
              Reaction = &
                   rr * &
                   Ions(iO_2DP_) * &
                   Neutrals(iN2_)

              IonSources(iN2P_)    = IonSources(iN2P_)   + Reaction
              NeutralSources(iO_3P_)  = NeutralSources(iO_3P_) + Reaction
              IonLosses(iO_2DP_)   = IonLosses(iO_2DP_)  + Reaction
              NeutralLosses(iN2_)  = NeutralLosses(iN2_) + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + Reaction * 0.859
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 0.491
              ChemicalHeatingS(iop2d_n2) =  &
                   ChemicalHeatingS(iop2d_n2) + &
                   Reaction * 1.33

              ! O+(2P) + N2 -> N2+ + O + 3.05 eV

              rr = 2.0e-16 * ti3m055
              Reaction = &
                   rr * &
                   Ions(iO_2PP_) * &
                   Neutrals(iN2_)

              IonSources(iN2P_)    = IonSources(iN2P_)   + Reaction
              NeutralSources(iO_3P_)  = NeutralSources(iO_3P_) + Reaction
              IonLosses(iO_2PP_)   = IonLosses(iO_2PP_)  + Reaction
              NeutralLosses(iN2_)  = NeutralLosses(iN2_) + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + Reaction * 1.941
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 1.109

               ChemicalHeatingS(iop2p_n2) =  &
                   ChemicalHeatingS(iop2p_n2) + &
                   Reaction * 3.02

              ! N2+ + O2 -> O2+ + N2 + 3.53 eV

               if (ti<=1000.0) then
                  rr = 5.1e-17 * ti3m116
               else 
                  rr = 1.26e-17 * ti10m067
               endif

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
                   ChemicalHeatingSub + Reaction * 1.822
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 1.708
              
              ChemicalHeatingS(in2p_o2) =  &
                   ChemicalHeatingS(in2p_o2) + &
                   Reaction * 3.53

              ! N2+ + O -> NO+ + N(2D) + 0.70 eV
              !!!!!             -> NO+ + N(4S) + 3.08 eV

              if (ti<=1500.0) then
                 rr = 1.33e-16 * ti3m044
              else 
                 rr = 6.55e-17 * ti15m023
              endif

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
                   Reaction * 0.477 

              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 0.223
                   
              ChemicalHeatingS(in2p_o) =  &
                   ChemicalHeatingS(in2p_o) + &
                   Reaction * 0.70
                   
              ! N2+ + e -> 2 N(2D) + 1.04 eV
              ! N2+ + e -> 2 N(4S) + 5.77 eV

              rr = 2.2e-13 * te3m039

              Reaction = &
                   rr * &
                   Ions(iN2P_) * &
                   Ions(ie_)

              NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + 2*Reaction*0.56
              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + 2*Reaction*0.44
              IonLosses(iN2P_)       = IonLosses(iN2P_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 1.04 *0.56 + &
                   Reaction * 5.77 * 0.44

              ChemicalHeatingS(in2p_e) =  &
                   ChemicalHeatingS(in2p_e) + &
                   Reaction * 3.44

              ! N2+ + N(4S) -> N2 + N+ + 2.48 eV

              rr = 1.0e-17

              Reaction = &
                   rr * &
                   Ions(iN2P_) * &
                   Neutrals(iN_4S_)

              NeutralSources(iN2_) = NeutralSources(iN2_) + Reaction
              IonSources(iNP_) = IonSources(iNP_) + Reaction
              NeutralLosses(iN_4S_) = NeutralLosses(iN_4S_) + Reaction
              IonLosses(iN2P_)       = IonLosses(iN2P_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 0.827
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 1.653
                  

              ! N2+ + O -> O+(4S) + N2 + 1.96 eV

              rr = 7.0e-18 * ti3m023

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
                   Reaction * 0.712
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 1.248

              ChemicalHeatingS(in2p_o) =  &
                   ChemicalHeatingS(in2p_o) + &
                   Reaction * 1.96

              ! N2+ + NO -> NO+ + N2 + 6.33 eV

              rr = 3.6e-16  ! Richards

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
                   Reaction * 3.274
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 3.056

              ! ----------------------------------------------------------
              ! O2+
              ! ----------------------------------------------------------

              ! -----------
              ! Solar EUV
              ! -----------

              Reaction = EuvIonRateS(iLon,iLat,iAlt,iO2P_,iBlock) * &
                         Neutrals(iO2_) + &
                         O2PERateS(iLon,iLat,iAlt,2,iBlock)*&
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

              if (ti<=900.0) then
                 rr = 1.6e-17 * ti3m052
              else
                 rr = 9e-18 * ti9m092
              endif

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
                   Reaction * 1.033
				  
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 0.517

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
                   Reaction * 3.243
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 1.622
              
              ! -----------
              ! O+(2P) + O2 -> O2+ + O + 6.54 eV
              ! -----------

              rr = 1.3e-16
              
              Reaction = &
                   rr * &
                   Ions(iO_2PP_) * &
                   Neutrals(iO2_)

              NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
              IonSources(iO2P_)   = IonSources(iO2P_)   + Reaction
              o2ptotal = o2ptotal + reaction
              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction
              IonLosses(iO_2PP_)  = IonLosses(iO_2PP_) + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 4.36
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 2.18

              ! -----------
              ! O+(2P) + O2 -> O+(4S) + O2 + 5.016 eV
              ! -----------

              rr = 1.3e-16

              Reaction = &
                   rr * &
                   Ions(iO_2PP_) * &
                   Neutrals(iO2_)

              NeutralSources(iO2_) = NeutralSources(iO2_) + Reaction
              IonSources(iO_4SP_)   = IonSources(iO_4SP_)   + Reaction
              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction
              IonLosses(iO_2PP_)  = IonLosses(iO_2PP_) + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 1.672
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 3.344

              ! -----------
              ! N+ + O2 -> O2+ + N(4S) + 2.5 eV
              ! -----------

              if (ti<=1000.0) then
                 rr = 1.925e-16 * ti3m045
              else 
                 rr = 3.325e-16
              endif

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
                  Reaction * 1.739
				  
             ChemicalHeatingSubI = &
                  ChemicalHeatingSubI + Reaction * 0.761
             
             ChemicalHeatingS(inp_o2) =  &
                  ChemicalHeatingS(inp_o2) + &
                  Reaction * 2.486
             
             ! -----------
             ! N+ + O2 -> O2+ + N(2D) + 0.1 eV
             ! -----------

             if (ti<=1000.0) then
                rr =  0.825e-16 * ti3m045
             else
                rr = 1.425e-16
             endif

             Reaction = &
                  rr * &
                  Ions(iNP_) * &
                  Neutrals(iO2_)

             NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + Reaction

             IonSources(iO2P_)      = IonSources(iO2P_)      + Reaction
             o2ptotal               = o2ptotal + reaction
             NeutralLosses(iO2_)    = NeutralLosses(iO2_)    + Reaction
             IonLosses(iNP_)        = IonLosses(iNP_)       + Reaction

             ChemicalHeatingSub = &
                  ChemicalHeatingSub + &
                  Reaction * 0.07
				  
             ChemicalHeatingSubI = &
                  ChemicalHeatingSubI + Reaction * 0.03

             ChemicalHeatingS(inp_o2) =  &
                  ChemicalHeatingS(inp_o2) + &
                  Reaction * 0.1

             ! -----------
             ! O2+ + e -> O(1D) + O(1D) + 3.06 eV
             ! O2+ + e -> O(3P) + O(1D) + 5.02 eV
             ! O2+ + e -> O(3P) + O(3P) + 6.99 eV
             ! -----------

             if (ti<=1200.0) then
                rr = 1.95e-13 * te3m07
             else
                rr = 7.39e-14 * te12m056
             endif

             Reaction = &
                  rr * &
                  Ions(iO2P_) * &
                  Ions(ie_)
           
             if (UseNeutralConstituent(iO_1D_)) then
                NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + 0.22*Reaction * 2.0
                NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + 0.42*Reaction
                NeutralSources(iO_1D_) = NeutralSources(iO_1D_) + 0.42*Reaction
                NeutralSources(iO_1D_) = NeutralSources(iO_1D_) + 0.31*Reaction * 2.0
                ! This really should be 0.05 to O(1D) and O(1S)
                NeutralSources(iO_1D_) = NeutralSources(iO_1D_) + 0.05*Reaction
                NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + 0.05*Reaction

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
                 
                NeutralSources(iO_3P_)    = &
                     NeutralSources(iO_3P_) + Reaction * 2.0

                ChemicalHeatingSub = &
                     ChemicalHeatingSub + &
                     Reaction * 5.0
                 
                ChemicalHeatingS(io2p_e) = &
                     ChemicalHeatingS(io2p_e) + &
                     Reaction * 5.0

             endif
              
             IonLosses(iO2P_)      = IonLosses(iO2P_)   + Reaction

             ! -----------
             ! O2+ + N(4S) -> NO+ + O + 4.21 eV
             ! -----------

             rr = 1.0e-16 ! Richards

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
                  Reaction * 2.746
				   
             ChemicalHeatingSubI = &
                  ChemicalHeatingSubI + Reaction * 1.464

             ChemicalHeatingS(io2p_n) =  &
                  ChemicalHeatingS(io2p_n) + &
                  Reaction * 4.25

             ! -----------
             ! O2+ + N(2D) -> NO+ + O + 6.519 eV
             ! -----------

             rr = 1.8e-16 ! Richards

             Reaction = &
                  rr * &
                  Ions(iO2P_) * &
                  Neutrals(iN_2D_)

             NeutralSources(iO_3P_)   = NeutralSources(iO_3P_)   + Reaction
             IonSources(iNOP_)     = IonSources(iNOP_)     + Reaction
             NeutralLosses(iN_2D_) = NeutralLosses(iN_2D_) + Reaction
             IonLosses(iO2P_)      = IonLosses(iO2P_)      + Reaction
             
             ChemicalHeatingSub = &
                  ChemicalHeatingSub + &
                  Reaction * 4.252
				   
             ChemicalHeatingSubI = &
                  ChemicalHeatingSubI + Reaction * 2.267

             ! -----------
             ! O2+ + N(2P) -> O2+ + N(4S) + 3.565 eV
             ! -----------

             rr = 2.2e-17 ! Richards

             Reaction = &
                  rr * &
                  Ions(iO2P_) * &
                  Neutrals(iN_2P_)

             NeutralSources(iN_4S_)   = NeutralSources(iN_4S_)   + Reaction
             IonSources(iO2P_)     = IonSources(iO2P_)     + Reaction
             NeutralLosses(iN_2P_) = NeutralLosses(iN_2P_) + Reaction
             IonLosses(iO2P_)      = IonLosses(iO2P_)      + Reaction
              
             ChemicalHeatingSub = &
                  ChemicalHeatingSub + &
                  Reaction * 2.479

             ChemicalHeatingSubI = &
                  ChemicalHeatingSubI + Reaction * 1.086

             ! -----------
             ! O2+ + NO -> NO+ + O2 + 2.813 eV
             ! -----------

              rr = 4.5e-16 ! schunk and nagy

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
                   Reaction * 1.361
				   
				   ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 1.452

              ChemicalHeatingS(io2p_no) =  &
                   ChemicalHeatingS(io2p_no) + &
                    Reaction * 2.813

!!! Temp change to stop crash
!
!              ! -----------
!              ! O2+ + N2 -> NO+ + NO + 0.9333 eV
!              ! -----------
!               
!              rr = 5.0e-22
!               
!              Reaction = &
!                   rr * &
!                   Ions(iO2P_) * &
!                   Neutrals(iN2_)
!
!              NeutralSources(iNO_) = NeutralSources(iNO_) + Reaction
!
!              IonSources(iNOP_)    = IonSources(iNOP_)    + Reaction
!              NeutralLosses(iN2_)  = NeutralLosses(iN2_)  + Reaction
!              IonLosses(iO2P_)     = IonLosses(iO2P_)     + Reaction
!
!              ChemicalHeatingSub = &
!                   ChemicalHeatingSub + &
!                   Reaction * 0.9333
!
!              ChemicalHeatingS(io2p_n2) =  &
!                   ChemicalHeatingS(io2p_n2) + &
!                   Reaction * 0.9333

              ! ----------------------------------------------------------
              ! O(4S)+
              ! ----------------------------------------------------------

              ! Solar EUV

!              Reaction = EuvIonRateS(iLon,iLat,iAlt,iO_4SP_,iBlock) * &
!                   Neutrals(iO_3P_)

              rr=EuvIonRateS(iLon,iLat,iAlt,iO_4SP_,iBlock) + &
                    OPERateS(iLon,iLat,iAlt,1,iBlock)
              Reaction = rr*Neutrals(iO_3P_)

              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction

              ! Aurora

              ! Aurora goes 0.4, 0.4, 0.2 into O(4S), O(2D) and O(2P) respectively
              Reaction = &
                   0.4 * AuroralIonRateS(iLon,iLat,iAlt,iO_3P_, iBlock) + &
                   0.4 * IonPrecipIonRateS(iLon,iLat,iAlt,iO_3P_, iBlock)

              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction

              ! Aurora goes 0.4, 0.4, 0.2 into O(3P), O(2D) and O(2P) respectively
              Reaction = &
                   0.4 * AuroralIonRateS(iLon,iLat,iAlt,iO_3P_, iBlock) + &
                   0.4 * IonPrecipIonRateS(iLon,iLat,iAlt,iO_3P_, iBlock)

              IonSources(iO_2DP_) = IonSources(iO_2DP_) + Reaction
              NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction

              ! Aurora goes 0.4, 0.4, 0.2 into O(3P), O(2D) and O(2P) respectively
              Reaction = &
                   0.2 * AuroralIonRateS(iLon,iLat,iAlt,iO_3P_, iBlock) + &
                   0.2 * IonPrecipIonRateS(iLon,iLat,iAlt,iO_3P_, iBlock)

              IonSources(iO_2PP_) = IonSources(iO_2PP_) + Reaction
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

              ! We create and loose the same amount of O 
              ! (when O(1D) is not used...
             
              if (UseNeutralConstituent(iO_1D_)) then
                 
                 NeutralSources(iO_3P_)  = &
                      NeutralSources(iO_3P_) + 0.5 * Reaction
                 NeutralSources(iO_1D_)  = &
                      NeutralSources(iO_1D_) + 0.5 * Reaction
                 
                 NeutralLosses(iO_3P_)   = NeutralSources(iO_3P_)  + Reaction

                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 1.655 * 0.5 + &
                      Reaction * 0.675 * 0.5
					  
                 ChemicalHeatingSubI = &
                      ChemicalHeatingSubI + &
                      Reaction * 1.655 * 0.5 + &
                      Reaction * 0.675 * 0.5 
              else

                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 1.655
					  
                 ChemicalHeatingSubI = &
                      ChemicalHeatingSubI + Reaction * 1.655

              endif
              
              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2DP_)  = IonLosses(iO_2DP_)  + Reaction
              
              ! -----------
              ! O+(2D) + e -> O+(4S) + e + 3.31 eV
              ! -----------

              rr = 6.03e-14 * te3m05

              Reaction = &
                   rr * &
                   Ions(iO_2DP_) * &
                   Ions(ie_)

              ! We create and loose the same amount of e
              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2DP_)  = IonLosses(iO_2DP_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 3.31
				   
              ChemicalHeatingS(iop2d_e) =  &
                   ChemicalHeatingS(iop2d_e) + &
                   Reaction * 0.0

!!! Temp change to stop crash
!               ! -----------
!               ! O+(2D) + N2 -> O+(4S) + N2 + 3.31 eV
!               ! -----------
!
!               rr = 8.0e-16
!
!               Reaction = &
!                    rr * &
!                    Ions(iO_2DP_) * &
!                    Neutrals(iN2_)
!
!               ! We create and loose the same amount of N2
!               IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
!               IonLosses(iO_2DP_)  = IonLosses(iO_2DP_)  + Reaction
!
!               ChemicalHeatingSubI = &
!                    ChemicalHeatingSubI + &
!                    Reaction * 3.31
!
!               ChemicalHeatingS(iop2d_n2) =  &
!                    ChemicalHeatingS(iop2d_n2) + &
!                    Reaction * 3.31

              ! -----------
              ! O+(2P) + O -> O+(4S) + O + 5.0 eV
              ! -----------

              rr = 4.0e-16

              Reaction = &
                   rr * &
                   Ions(iO_2PP_) * &
                   Neutrals(iO_3P_)

              ! We create and loose the same amount of O
              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2PP_)  = IonLosses(iO_2PP_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 2.5
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 2.5

              ChemicalHeatingS(iop2p_o) =  &
                   ChemicalHeatingS(iop2p_o) + &
                   Reaction * 5.0

              ! -----------
              ! O+(2P) + e -> O+(4S) + e + 5.0 eV
              ! -----------

              rr = 3.03e-14 * te3m05

              Reaction = &
                   rr * &
                   Ions(iO_2PP_) * &
                   Ions(ie_)

              ! We create and loose the same amount of e
              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2PP_)  = IonLosses(iO_2PP_)  + Reaction

              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 5.0

              ! maybe should be 0 because the energy should go the electrons
              ChemicalHeatingS(iop2p_e) =  &
                   ChemicalHeatingS(iop2p_e) + &
                   Reaction*5 

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

              ! -----------
              ! N+ + O2 -> O+(4S) + NO + 2.31 eV
              ! -----------
              
              if (ti<=1000.0) then
                 rr = 0.275e-16 * ti3m045
              else 
                 rr = 0.475e-16
              endif
			 
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
                   Reaction * 0.803
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 1.507

              ChemicalHeatingS(inp_o2) =  &
                   ChemicalHeatingS(inp_o2) + &
                   Reaction * 2.31

              ! -----------
              ! O+(4S) + N2 -> NO+ + N(4S) + 1.10 eV
              ! -----------

              if (ti<=1000.0) then
                 rr = 1.2e-18 * ti3m045
              else
                 rr = 7.0e-19 * ti10m212
              endif

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
                   Reaction * 0.75

              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 0.35

              ChemicalHeatingS(iop_n2) =  &
                   ChemicalHeatingS(iop_n2) + &
                   Reaction * 1.10

              ! -----------
              ! O+(4S) + NO -> NO+ + O + 4.36 eV
              ! -----------

              rr = 7.0e-19 * ti3m087

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
                   Reaction * 2.844
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 1.516

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
                   Reaction * 0.677
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 0.773
              
              ! ----------------------------------------------------------
              ! O(2D)+
              ! ----------------------------------------------------------

              ! Solar EUV

!              Reaction = EuvIonRateS(iLon,iLat,iAlt,iO_2DP_,iBlock) * &
!                   Neutrals(iO_3P_)

              rr=EuvIonRateS(iLon,iLat,iAlt,iO_2DP_,iBlock) + &
                    OPERateS(iLon,iLat,iAlt,2,iBlock)
              Reaction = rr * Neutrals(iO_3P_)

              IonSources(iO_2DP_) = IonSources(iO_2DP_) + Reaction
              NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction

              ! -----------
              ! O+(2P) + e -> O+(2D) + e + 1.69 eV
              ! -----------

              rr = 1.84e-13 * te3m05

              Reaction = &
                   rr * &
                   Ions(iO_2PP_) * &
                   Ions(ie_)

              ! We create and loose the same amount of e
              IonSources(iO_2DP_) = IonSources(iO_2DP_) + Reaction
              IonLosses(iO_2PP_)  = IonLosses(iO_2PP_)  + Reaction

              ChemicalHeatingSubI = &
                   ChemicalHeatingSub + &
                   Reaction * 1.69
				   
              ChemicalHeatingS(iop2p_e) =  &
                   ChemicalHeatingS(iop2p_e) + &
                   Reaction * 1.69
                   
               ! -----------
              ! O+(2D) + N2 -> NO+ + N + 4.41 eV
              ! -----------

              rr = 2.5e-17

              Reaction = &
                   rr * &
                   Ions(iO_2DP_) * &
                   Neutrals(iN2_)

              NeutralLosses(iN2_)  = NeutralLosses(iN2_)  + Reaction
              IonSources(iNOP_) = IonSources(iNOP_) + Reaction
              IonLosses(iO_2DP_)  = IonLosses(iO_2DP_)  + Reaction
              NeutralSources(iN_4S_)  = NeutralSources(iN_4S_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 3.007
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 1.403
                   
              ! -----------
              ! O+(2D) + NO -> NO+ + O + 4.37 eV
              ! -----------

              rr = 1.2e-15

              Reaction = &
                   rr * &
                   Ions(iO_2DP_) * &
                   Neutrals(iNO_)

              NeutralLosses(iNO_)  = NeutralLosses(iNO_)  + Reaction
              IonSources(iNOP_) = IonSources(iNOP_) + Reaction
              IonLosses(iO_2DP_)  = IonLosses(iO_2DP_)  + Reaction
              NeutralSources(iO_3P_)  = NeutralSources(iO_3P_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 2.85
				  
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 1.52
              
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

              !Reaction = EuvIonRateS(iLon,iLat,iAlt,iO_2PP_,iBlock) * &
              !     Neutrals(iO_3P_)
!              rr=EuvIonRateS(iLon,iLat,iAlt,iO_2PP_,iBlock) 

              rr=EuvIonRateS(iLon,iLat,iAlt,iO_2PP_,iBlock) + &
                   OPERateS(iLon,iLat,iAlt,3,iBlock)
              Reaction = rr*Neutrals(iO_3P_)

              IonSources(iO_2PP_) = IonSources(iO_2PP_) + Reaction
              NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction

              ! ----------------------------
              ! Atomic N Photoionization
              ! ----------------------------
              ! ----------------------------------------------------------
              ! N(4S) + hv -> N+
              ! ----------------------------------------------------------
              Reaction = EuvIonRateS(iLon,iLat,iAlt,iNP_,iBlock) * &
                   Neutrals(iN_4S_)

              IonSources(iNP_) = IonSources(iNP_) + Reaction
              NeutralLosses(iN_4S_)  = NeutralLosses(iN_4S_)  + Reaction

!              ! ----------------------------
!              ! Atomic He Photoionization
!              ! ----------------------------
!              ! ----------------------------------------------------------
!              ! He + hv --> He+  + e-
!              ! ----------------------------------------------------------
!              Reaction = EuvIonRateS(iLon,iLat,iAlt,iHeP_,iBlock) * &
!                   Neutrals(iHe_)
!
!              NeutralLosses(iHe_) = NeutralLosses(iHe_) + Reaction
!              IonSources(iHeP_) = IonSources(iHeP_) + Reaction

! ----------------------------
! NO Photoionization

!              IonSources(iO_2PP_) = IonSources(iO_2PP_) + Reaction
!              NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction

!!! Temp change to stop crash
!               ! -----------
!               ! O+(2P) + N2 -> N+ + NO + 0.70 eV
!               ! -----------
!
!               rr = 1.0e-16
!
!               Reaction = &
!                    rr * &
!                    Ions(iO_2PP_) * &
!                    Neutrals(iN2_)
!
!               NeutralSources(iNO_) = NeutralSources(iNO_) + Reaction
!               IonSources(iNP_)     = IonSources(iNP_)     + Reaction
!               NeutralLosses(iN2_)  = NeutralLosses(iN2_)  + Reaction
!               IonLosses(iO_2PP_)   = IonLosses(iO_2PP_)   + Reaction
!
!               ChemicalHeatingSub = &
!                    ChemicalHeatingSub + &
!                    Reaction * 0.70
!
!               ChemicalHeatingS(iop2p_n2) =  &
!                    ChemicalHeatingS(iop2p_n2) + &
!                    Reaction * 0.70

               ! ----------------------------------------------------------
               ! N+
               ! ----------------------------------------------------------
               
              !! Temp change to stop crash
              ! ! -----------
              ! ! O2+ + N(2D) -> N+ + O2 + 0.0 eV
              ! ! -----------
              !
              ! rr = 8.65e-17
              !
              ! Reaction = &
              !      rr * &
              !      Ions(iO2P_) * &
              !      Neutrals(iN_2D_)
              !
              ! NeutralSources(iO2_)  = NeutralSources(iO2_)  + Reaction
              ! IonSources(iNP_)      = IonSources(iNP_)      + Reaction
              ! NeutralLosses(iN_2D_) = NeutralLosses(iN_2D_) + Reaction
              ! IonLosses(iO2P_)      = IonLosses(iO2P_)      + Reaction
              ! ChemicalHeatingSub = &
              !      ChemicalHeatingSub + &
              !      Reaction * 0.0

!              ! -----------
!              ! Shunk and Nagy R29
!              ! He+ + N2 -> N+ + N + He + 0.28 eV
!              ! -----------
!
!              rr = 7.8e-10/1.0e6
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
!
!              ! -----------
!              ! Shunk and Nagy R30
!              ! He+ + N2 -> N2+ + He ( + ??? eV)
!              ! -----------
!
!              rr = 5.2e-10/1.0e6
!
!              Reaction = &
!                   rr * &
!                   Ions(iHeP_) * &
!                   Neutrals(iN2_)
!
!              NeutralSources(iHe_)   = NeutralSources(iHe_)   + Reaction
!              IonSources(iN2P_)      = IonSources(iN2P_)      + Reaction
!              NeutralLosses(iN2_)    = NeutralLosses(iN2_)    + Reaction
!              IonLosses(iHeP_)       = IonLosses(iHeP_)       + Reaction
!
!              !ChemicalHeatingSub = &
!              !     ChemicalHeatingSub + &
!              !     Reaction * 0.28
!
!              ! -----------
!              ! Shunk and Nagy R31
!              ! He+ + O2 -> O+ O + He ( + ??? eV)
!              ! -----------
!
!              rr = 9.7e-10/1.0e6
!
!              Reaction = &
!                   rr * &
!                   Ions(iHeP_) * &
!                   Neutrals(iO2_)
!
!              NeutralSources(iHe_) = NeutralSources(iHe_) + Reaction
!              NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
!              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
!              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction
!              IonLosses(iHeP_) = IonLosses(iHeP_) + Reaction
!
!              ! -----------
!              ! Shunk and Nagy Radiative Recombination
!              ! He+ + e- -> He ( + ??? eV)
!              ! -----------
!
!              rr = 4.8e-12/1.0e6 * te07
!
!              Reaction = &
!                   rr * &
!                   Ions(iHeP_) * &
!                   Ions(ie_)
!
!              NeutralSources(iHe_) = NeutralSources(iHe_) + Reaction
!              IonLosses(iHeP_) = IonLosses(iHeP_) + Reaction

              !ChemicalHeatingSub = &
              !     ChemicalHeatingSub + &
              !     Reaction * 0.28

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
              ! N+ + NO --> N2+ + O + 2.2 eV
              ! -----------

               rr = 8.33e-17 * ti3m024

               Reaction = &
                    rr * &
                    Ions(iNP_) * &
                    Neutrals(iNO_)

               NeutralSources(iO_3P_)   = NeutralSources(iO_3P_)   + Reaction
               IonSources(iN2P_)      = IonSources(iN2P_)      + Reaction
               NeutralLosses(iNO_) = NeutralLosses(iNO_) + Reaction
               IonLosses(iNP_)    = IonLosses(iNP_)    + Reaction

               ChemicalHeatingSub = &
                    ChemicalHeatingSub + &
                    Reaction * 1.4
				   
               ChemicalHeatingSubI = &
                    ChemicalHeatingSubI + Reaction * 0.8

               ! -----------
               ! N+ + NO --> NO+ + N(4S) + 3.4 eV
               ! -----------

               rr = 4.72e-16 * ti3m024

               Reaction = &
                    rr * &
                    Ions(iNP_) * &
                    Neutrals(iNO_)

               NeutralSources(iN_4S_)   = NeutralSources(iN_4S_)   + Reaction
               IonSources(iNOP_)      = IonSources(iNOP_)      + Reaction
               NeutralLosses(iNO_) = NeutralLosses(iNO_) + Reaction
               IonLosses(iNP_)    = IonLosses(iNP_)    + Reaction

               ChemicalHeatingSub = &
                    ChemicalHeatingSub + &
                    Reaction * 2.318
				   
               ChemicalHeatingSubI = &
                    ChemicalHeatingSubI + Reaction * 1.082

               ! -----------
               ! N+ + O2 --> NO+ + O(3P) + 6.67 eV
               ! -----------
              
               if (ti<=1000) then
                  rr = 0.495e-16 * ti3m045
               else 
                  rr = 0.855e-16
               endif

               Reaction = &
                    rr * &
                    Ions(iNP_) * &
                    Neutrals(iO2_)

               NeutralSources(iO_3P_)   = NeutralSources(iO_3P_)   + Reaction
               IonSources(iNOP_)      = IonSources(iNOP_)      + Reaction
               NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction
               IonLosses(iNP_)    = IonLosses(iNP_)    + Reaction

               ChemicalHeatingSub = &
                    ChemicalHeatingSub + &
                    Reaction * 4.35
				   
               ChemicalHeatingSubI = &
                    ChemicalHeatingSubI + Reaction * 2.32

              ! -----------
              ! O+(2D) + N -> N+ + O + 1.0 eV
              ! -----------

              rr = 1.5e-16

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
                   Reaction * 0.467
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 0.533

              ! -----------
              ! N+ + O2 -> NO+ + O(1D) + 4.71 eV
              ! -----------

              if (ti<=1000.0) then
                 rr = 1.98e-16 * ti3m045
              else 
                 rr = 3.42e-16
              endif

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
                   Reaction * 3.072
              
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 1.638

              ChemicalHeatingS(inp_o2) =  &
                   ChemicalHeatingS(inp_o2) + &
                   Reaction * 6.67

              ! -----------
              ! N+ + O -> O+ + N + 0.93 eV
              ! -----------

              rr = 2.2e-18

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
                   Reaction * 0.496
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 0.434

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
              ! NO+ + e -> O + N(2D) + 0.38 eV 
              ! -----------

              rr = 3.4e-13 * te3m085

              Reaction = &
                   rr * &
                   Ions(iNOP_) * &
                   Ions(ie_)

              NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + Reaction
              NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
              IonLosses(iNOP_)       = IonLosses(iNOP_) + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 0.38

              ChemicalHeatingS(inop_e) =  &
                   ChemicalHeatingS(inop_e) + &
                   Reaction * 0.38
                   
              ! -----------
              ! NO+ + e -> O + N(4S) + 2.77 eV
              ! -----------

              rr = 0.6e-13 * te3m085 

              Reaction = &
                   rr * &
                   Ions(iNOP_) * &
                   Ions(ie_)


              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
              IonLosses(iNOP_)       = IonLosses(iNOP_)       + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 2.77

               ChemicalHeatingS(inop_e) =  &
                   ChemicalHeatingS(inop_e) + &
                   Reaction * 2.77

              ! ----------------------------------------------------------
              ! N(4S)
              ! ----------------------------------------------------------

              ! -----------
              ! N(2D) + e -> N(4S) + e + 2.38 eV
              ! -----------

              rr = 3.86e-16 * te3m081

              Reaction = &
                   rr * &
                   Neutrals(iN_2D_) * &
                   Ions(ie_)

              ! We create and loose the same amount of e
              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              NeutralLosses(iN_2D_)  = NeutralLosses(iN_2D_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 2.38
				   
              ! -----------
              ! N(2D) + O -> N(4S) + O(3P) + 2.38 eV
              ! N(2D) + O -> N(4S) + O(1D) + 0.42 eV
              ! -----------

              rr = 6.9e-19

              Reaction = &
                   rr * &
                   Neutrals(iN_2D_) * &
                   Neutrals(iO_3P_)

              if (UseNeutralConstituent(iO_1D_)) then
                 NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + 0.9*Reaction
                 NeutralSources(iO_1D_) = NeutralSources(iO_1D_) + 0.1*Reaction
                 NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction
                 
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

              rr=4.5e-6*exp(-1.e-8*(Neutrals(iO2_)*1.e-6)**0.38)

              Reaction = &
                   rr * &
                   Neutrals(iNO_)

              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
              NeutralLosses(iNO_)    = NeutralLosses(iNO_)    + Reaction

              ! -----------
              ! N(4S) + O2 -> NO + O + 1.385 eV
              ! -----------
              
              rr = 1.5e-20 * tn * exp(-3270/tn) 

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
              rr = 3.4e-17


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

              rr = 9.7e-18 * exp(-185/tn)

              Reaction = &
                   rr * &
                   Neutrals(iN_2D_) * &
                   Neutrals(iO2_)

              if (UseNeutralConstituent(iO_1D_)) then
                 
                 NeutralSources(iO_3P_) = &
                      NeutralSources(iO_3P_) + 0.9 * Reaction
                 NeutralSources(iO_1D_) = &
                      NeutralSources(iO_1D_) + 0.1 * Reaction
                 
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

              rr = 6.7e-17

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

!!! Temp change to stop crash
!              ! ----------------------------------------------------------
!              ! N(2P)
!              ! ----------------------------------------------------------
!			  
!              ! -----------
!              ! N(2P) + e -> N(2D) + e + 1.19 eV
!              ! -----------
!
!              rr = 9.5e-15
!
!              Reaction = &
!                   rr * &
!                   Neutrals(iN_2P_) * &
!                   Ions(ie_)
!
!              ! We create and loose the same amount of e
!              NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + Reaction
!              NeutralLosses(iN_2P_)  = NeutralLosses(iN_2P_)  + Reaction
!
!              ChemicalHeatingSub = &
!                   ChemicalHeatingSub + &
!                   Reaction * 1.19
!				   
!              ! -----------
!              ! N(2P) + e -> N(4S) + e + 3.57 eV
!              ! -----------
!
!              rr = 2.04e-16 * (te3m085**(-1))
!
!              Reaction = &
!                   rr * &
!                   Neutrals(iN_2P_) * &
!                   Ions(ie_)
!
!              ! We create and loose the same amount of e
!              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
!              NeutralLosses(iN_2P_)  = NeutralLosses(iN_2P_)  + Reaction
!
!              ChemicalHeatingSub = &
!                   ChemicalHeatingSub + &
!                   Reaction * 3.57
!				   
!              ! -----------
!              ! N(2P) + NO -> N(4S) + NO + 3.44 eV
!              ! -----------
!
!              rr = 1.8e-16
!
!              Reaction = &
!                   rr * &
!                   Neutrals(iN_2P_) * &
!                   Neutrals(iNO_)
!
!              NeutralSources(iN_4S_)  = NeutralSources(iN_4S_)  + Reaction
!              NeutralLosses(iN_2P_) = NeutralLosses(iN_2P_) + Reaction
!
!              ChemicalHeatingSub = &
!                   ChemicalHeatingSub + &
!                   Reaction * 3.44
!				   
!              ! -----------
!              ! N(2P) + O(3P) -> N(2D) + O(3P) + 1.19 eV
!              ! -----------
!
!              rr = 1.7e-17
!
!              Reaction = &
!                   rr * &
!                   Neutrals(iN_2P_) * &
!                   Neutrals(iO_3P_)
!
!              NeutralSources(iN_2D_)  = NeutralSources(iN_2D_)  + Reaction
!              NeutralLosses(iN_2P_) = NeutralLosses(iN_2P_) + Reaction
!
!              ChemicalHeatingSub = &
!                   ChemicalHeatingSub + &
!                   Reaction * 1.19
!				   
!              ! -----------
!              ! N(2P) + O2 -> NO + O(3P) + 4.95 eV
!              ! -----------
!
!              rr = 3.09e-18 * exp(-60/Tn)
!
!              Reaction = &
!                   rr * &
!                   Neutrals(iN_2P_) * &
!                   Neutrals(iO2_)
!
!              NeutralSources(iNO_)  = NeutralSources(iNO_)  + Reaction
!              NeutralSources(iO_3P_)  = NeutralSources(iO_3P_)  + Reaction
!              NeutralLosses(iN_2P_) = NeutralLosses(iN_2P_) + Reaction
!              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction
!
!              ChemicalHeatingSub = &
!                   ChemicalHeatingSub + &
!                   Reaction * 4.95

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
                 ! O(1D) + e -> O(3P) + e + 1.96 eV
                 ! ------------                                             
                 
                 rr = 2.87e-16 * te3m091
                 
                 Reaction = &
                      rr * &
                      Neutrals(iO_1D_) * Ions(ie_)
                 
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
                 
                 rr = 1.8e-17 * exp(107/Tn)

                 Reaction = &
                      rr * &
                      Neutrals(iO_1D_) * &
                      Neutrals(iN2_)
                 
                 ! We create and loose the same amount of N2
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
                 
                 rr = 3.2e-17 * exp(67/Tn)
                 
                 Reaction = &
                      rr * &
                      Neutrals(iO_1D_) * &
                      Neutrals(iO2_)
                 
                 ! We create and loose the same amount of O2
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
                 
                 rr = 6.47e-18 * ((Tn/300)**0.14)
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
                   (Neutrals(iO2_)*1.e-6)**0.8855)*szap
              Reaction = &
                   rr * &
                   Neutrals(iNO_)

              IonSources(iNOP_)   = IonSources(iNOP_)   + Reaction
              NeutralLosses(iNO_) = NeutralLosses(iNO_) + Reaction
              
              !---- Ions

              if (.not. UseIonChemistry) then
                 IonSources = 0.0
                 IonLosses = 0.0
              else
                 do iIon = 1, nIons-1
                    if (.not.UseIonConstituent(iIon)) then
                       IonSources(iIon) = 0.0
                       IonLosses(iIon) = 0.0
                    endif
                 enddo
              endif

              if (.not. UseNeutralChemistry) then
                 NeutralSources = 0.0
                 NeutralLosses = 0.0
              else
                 do iNeutral = 1, nSpeciesTotal
                    if (.not.UseNeutralConstituent(iNeutral)) then
                       NeutralSources(iNeutral) = 0.0
                       NeutralLosses(iNeutral) = 0.0
                    endif
                 enddo
              endif
              
              ! Take Implicit time step
              Ions(ie_) = 0.0
              do iIon = 1, nIons-1
                 ionso = IonSources(iIon)
                 ionlo = IonLosses(iIon)/(Ions(iIon)+1.0e-6)
                 Ions(iIon) = (Ions(iIon) + ionso * DtSub) / &
                      (1 + DtSub * ionlo)
                 ! sum for e-
                 Ions(ie_) = Ions(ie_) + Ions(iIon)
              enddo

              do iNeutral = 1, nSpeciesTotal

                 neuso = NeutralSources(iNeutral)
                 neulo = NeutralLosses(iNeutral) / Neutrals(iNeutral)

                 Neutrals(iNeutral)=(Neutrals(iNeutral) + neuso * DtSub) / &
                      (1 + DtSub * neulo)

                 NeutralSourcesTotal(ialt,iNeutral) = &
                      NeutralSourcesTotal(ialt,iNeutral) + &
                      NeutralSources(iNeutral) * DtSub

                 NeutralLossesTotal(ialt,iNeutral) = &
                      NeutralLossesTotal(ialt,iNeutral) + &
                      NeutralLosses(iNeutral) * DtSub
                 
              enddo

              ChemicalHeatingRate(iLon,iLat,iAlt) = &
                   ChemicalHeatingRate(iLon,iLat,iAlt) + &
                   ChemicalHeatingSub * DtSub
				   
              ChemicalHeatingRateIon(iLon,iLat,iAlt) = &
                   ChemicalHeatingRateIon(iLon,iLat,iAlt) + &
                   ChemicalHeatingSubI * DtSub

              ChemicalHeatingRateEle(iLon,iLat,iAlt) = &
                   ChemicalHeatingRateEle(iLon,iLat,iAlt) + &
                   ChemicalHeatingSubE * DtSub

              ChemicalHeatingSpecies(iLon,iLat,iAlt,:) = &
                   ChemicalHeatingSpecies(iLon,iLat,iAlt,:) + &
                   ChemicalHeatingS * DtSub

              EmissionTotal = EmissionTotal + Emission(:)*DtSub
              
              DtTotal = DtTotal + DtSub

              if (DtSub < DtMin) DtMin = DtSub

              if (DtSub < 1.0e-9 .and. abs(DtTotal-Dt) > DtSub) then
                 write(*,*) "Chemistry is too fast!!", DtSub

                 ! Check Ions
                 do iIon = 1, nIons
                    write(*,*) "Ion Source/Loss : ", &
                         iIon, IonSources(iIon), IonLosses(iIon)
                 enddo
                 do iNeutral = 1, nSpeciesTotal
                    write(*,*) "Neutral Source/Loss : ", iAlt, &
                         iNeutral, NeutralSources(iNeutral), &
                         NeutralLosses(iNeutral), Neutrals(iNeutral)
                 enddo

                 call stop_gitm("Chemistry is too fast!!")
              endif

              nIters = nIters + 1

           enddo

           IDensityS(iLon,iLat,iAlt,:,iBlock) = Ions
           NDensityS(iLon,iLat,iAlt,:,iBlock) = Neutrals

           Emissions(iLon, iLat, iAlt, :, iBlock) =  &
                Emissions(iLon, iLat, iAlt, :, iBlock) + EmissionTotal

        enddo
     enddo
  enddo
 
  if (iDebugLevel > 3) then
     do iIon = 1, nIons
        write(*,*) "====> calc_chemistry: Max Ion Density: ", iIon, &
             maxval(IDensityS(1:nLons,1:nLats,(nAlts*4)/5,iIon,iBlock))
     enddo
  endif

  if (iDebugLevel > 2) &
       write(*,*) "===> calc_chemistry: Average Dt for this timestep : ", &
       (Dt*nLats*nLons*nAlts)/nIters

  call end_timing("calc_chemistry")

end subroutine calc_chemistry
