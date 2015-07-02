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
  real :: AerosolLosses(nSpeciesTotal)
  real :: H2AerosolProd
  real :: IsoScalingFactor
  real :: DtSub, DtOld, DtTotal, DtMin, DtAve, Source, Reaction, tr, tr3, rr
  real :: te3, ti, te, tn, tn1, tn06, dtsubtmp, losstmp, dentmp, l, t,m1,m2,y1,y2,k1,k2
  real :: Ions(nIons), Neutrals(nSpeciesTotal)
  real :: tli(nIons), tsi(nIons), tln(nSpeciesTotal), tsn(nSpeciesTotal)

  integer :: iLon, iLat, iAlt, iIon, nIters, iNeutral
  
  real :: lon, lat, alt
  real, dimension(1:2) :: msis_temp
  real, dimension(1:8) :: msis_dens
  real :: LonDeg, LatDeg, AltKm, LST
  real,dimension(7)    :: AP  

  real :: ChemicalHeatingSub,percent,o2total,o2ptotal
  real :: Emission(nEmissions), EmissionTotal(nEmissions)
  real, dimension(nLons,nLats,nAlts) :: &
       tr3d, tr33d,te12d, tr3m0443d, tr3m083d, tr3m043d, &
       te3m0393d, te3m0853d, te33d, te3m053d,te3m073d,te12m0563d,te227d
  
  real :: te3m05,te3m07,te12m056, tr3m044, tr3m04, tr3m08, te3m039, te3m085, rr_opn2, &
       te22m05

  real :: rt300te
  real :: rt300te03, rt300te07

  real :: rt3CH2_H
  real :: rtCH_CH4
  real :: rtC2H4_CH

  real :: rtCH3_CH30, rtCH3_CH3inf

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

  DtMin = Dt

  if (.not.UseIonChemistry) return

  call report("Chemistry",2)
  call start_timing("calc_chemistry")

  DtAve = 0.0

  nIters=0

  tr3d = (iTemperature(1:nLons,1:nLats,1:nAlts,iBlock) &
         + Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*&
         TempUnit(1:nLons,1:nLats,1:nAlts)) / 2.0


  do iLon = 1, nLons
     do iLat = 1, nLats
        do iAlt = 1, nAlts

           NeutralSourcesTotal = 0.0
           NeutralLossesTotal = 0.0

           te = eTemperature(iLon,iLat,iAlt,iBlock)
           ti = iTemperature(iLon,iLat,iAlt,iBlock)
           tn = Temperature(iLon,iLat,iAlt,iBlock)*&
                TempUnit(iLon,iLat,iAlt)

!!
!! Reaction Rates with Electron Temperature
          rt300te   = sqrt(300.0/te)
          rt300te07 = (300.0/te)**0.70
          rt300te03 = (300.0/te)**0.30
!!
!! Reaction Rates with Neutral Temperature
          rt3CH2_H   = exp(-370.0/tn)
          rtCH_CH4   = (tn**(-1.04))*exp(-36.1/tn)
          rtC2H4_CH  = (tn**(-0.546))*exp(-29.6/tn)

           DtTotal = 0.0
           EmissionTotal = 0.0

           Ions = IDensityS(iLon,iLat,iAlt,:,iBlock)

           Neutrals = NDensityS(iLon,iLat,iAlt,:,iBlock)
           AerosolLosses = AerosolTrappingLoss(iLon,iLat,iAlt,:,iBlock)
!           H2AerosolProd = H2AerosolProducgion(iLon,iLat,iAlt,iBlock)
           IsoScalingFactor = IsotopeScaling(iLon,iLat,iAlt,iBlock)

           niters = 0

           do while (DtTotal < Dt)
              
              ChemicalHeatingSub = 0.0
              ChemicalHeatingS = 0
              Emission = 0.0

              DtSub = Dt - DtTotal

              IonSources = 0.0
              NeutralSources = 0.0
              IonLosses  = 0.0
              NeutralLosses = 0.0

              !--------------------------------------------->>>>>
!\
! Nitrogen Photochemistry:--------------------------------------------------+
!/
              ! ----------------------------------------------------------
              ! N2 + hv ==> N(4S) + N(2D) 
              ! ----------------------------------------------------------
              !rr=EuvDissRateS(iLon,iLat,iAlt,iN4S_,iBlock) 
              rr=EuvDissRateS(iLon,iLat,iAlt,iN4S_,iBlock) + MagDissRateS(iLon,iLat,iAlt,iN4S_,iBlock)

              Reaction = rr * Neutrals(iN2_)
              ! This gives the chemistry in terms of kg/m^3/s

              NeutralLosses(iN2_) = NeutralLosses(iN2_) + Reaction
              NeutralSources(iN4S_) = NeutralSources(iN4S_) + Reaction

              ! 15N-14N
              rr=EuvDissRateS(iLon,iLat,iAlt,iN4S_,iBlock)*IsoScalingFactor + MagDissRateS(iLon,iLat,iAlt,iN4S_,iBlock)
              Reaction = rr * Neutrals(i15N2_)
              NeutralLosses(i15N2_) = NeutralLosses(i15N2_) + Reaction

              ! ----------------------------------------------------------
              ! N2 + hv ==> N(4S) + N+ 
              ! ----------------------------------------------------------
              ! rr = EuvIonRateS(iLon,iLat,iAlt,iNP_,iBlock)
              rr = EuvIonRateS(iLon,iLat,iAlt,iNP_,iBlock) + MagIonRateS(iLon,iLat,iAlt,iNP_,iBlock)

              Reaction = rr * Neutrals(iN2_)

              NeutralLosses(iN2_) = NeutralLosses(iN2_) + Reaction
              NeutralSources(iN4S_) = NeutralSources(iN4S_) + Reaction
              IonSources(iNP_) = IonSources(iNP_) + Reaction

              ! Isotope Losses

              rr=EuvIonRateS(iLon,iLat,iAlt,iNP_,iBlock)*IsoScalingFactor + MagIonRateS(iLon,iLat,iAlt,iNP_,iBlock)
              Reaction = rr * Neutrals(i15N2_)
              NeutralLosses(i15N2_) = NeutralLosses(i15N2_) + Reaction

              ! ----------------------------------------------------------
              ! N2 + hv ==> N2+ 
              ! ----------------------------------------------------------
              ! rr=EuvIonRateS(iLon,iLat,iAlt,iN2P_,iBlock)
              rr=EuvIonRateS(iLon,iLat,iAlt,iN2P_,iBlock) + MagIonRateS(iLon,iLat,iAlt,iN2P_,iBlock)

              Reaction = rr * Neutrals(iN2_)

              NeutralLosses(iN2_) = NeutralLosses(iN2_) + Reaction
              IonSources(iN2P_) = IonSources(iN2P_) + Reaction

              ! Isotope Losses
              rr=EuvIonRateS(iLon,iLat,iAlt,iN2P_,iBlock)*IsoScalingFactor + MagIonRateS(iLon,iLat,iAlt,iN2P_,iBlock)
              Reaction = rr * Neutrals(i15N2_)
              NeutralLosses(i15N2_) = NeutralLosses(i15N2_) + Reaction

!\
! Methane Photochemistry:-----------------------------------------------------+
!/
              ! ----------------------------------------------------------
              ! CH4 + hv ==> CH + H2 + H 
              ! ----------------------------------------------------------
              !rr=EuvDissRateS(iLon,iLat,iAlt,iCH_,iBlock)
              rr=EuvDissRateS(iLon,iLat,iAlt,iCH_,iBlock)+  MagDissRateS(iLon,iLat,iAlt,iCH_,iBlock)

              Reaction = rr * Neutrals(iCH4_)

              NeutralLosses(iCH4_) = NeutralLosses(iCH4_) + Reaction

              NeutralSources(iCH_) = NeutralSources(iCH_) + Reaction
              NeutralSources(iH2_) = NeutralSources(iH2_) + Reaction
              NeutralSources(iH_) = NeutralSources(iH_) + Reaction

              ! ----------------------------------------------------------
              ! CH4 + hv ==> 3CH2 + 2H 
              ! ----------------------------------------------------------
              !rr=EuvDissRateS(iLon,iLat,iAlt,i3CH2_,iBlock)

              rr=EuvDissRateS(iLon,iLat,iAlt,i3CH2_,iBlock) +  MagDissRateS(iLon,iLat,iAlt,i3CH2_,iBlock)

              Reaction = rr * Neutrals(iCH4_)

              NeutralLosses(iCH4_) = NeutralLosses(iCH4_) + Reaction

              NeutralSources(i3CH2_) = NeutralSources(i3CH2_) + Reaction
              NeutralSources(iH_) = NeutralSources(iH_) + 2*Reaction

              ! ----------------------------------------------------------
              ! CH4 + hv ==> 1CH2 + H2 
              ! ----------------------------------------------------------
              !rr=EuvDissRateS(iLon,iLat,iAlt,i1CH2_,iBlock)

              rr=EuvDissRateS(iLon,iLat,iAlt,i1CH2_,iBlock) +  MagDissRateS(iLon,iLat,iAlt,i1CH2_,iBlock)

              Reaction = rr * Neutrals(iCH4_)

              NeutralLosses(iCH4_) = NeutralLosses(iCH4_) + Reaction

              NeutralSources(i1CH2_) = NeutralSources(i1CH2_) + Reaction
              NeutralSources(iH2_) = NeutralSources(iH2_) + Reaction

              ! ----------------------------------------------------------
              ! CH4 + hv ==> CH3 + H 
              ! ----------------------------------------------------------
              ! rr=EuvDissRateS(iLon,iLat,iAlt,iCH3_,iBlock)
              rr=EuvDissRateS(iLon,iLat,iAlt,iCH3_,iBlock) +  MagDissRateS(iLon,iLat,iAlt,iCH3_,iBlock)
              Reaction = rr * Neutrals(iCH4_) 

              NeutralLosses(iCH4_) = NeutralLosses(iCH4_) + Reaction

              NeutralSources(iCH3_) = NeutralSources(iCH3_) + Reaction
              NeutralSources(iH_) = NeutralSources(iH_) + Reaction

              ! ----------------------------------------------------------
              ! CH4 + hv ==> CH3+ + H 
              ! ----------------------------------------------------------
              ! rr=EuvIonRateS(iLon,iLat,iAlt,iCH3P_,iBlock)
              rr=EuvIonRateS(iLon,iLat,iAlt,iCH3P_,iBlock) + MagIonRateS(iLon,iLat,iAlt,iCH3P_,iBlock)

              Reaction = rr * Neutrals(iCH4_) 

              NeutralLosses(iCH4_) = NeutralLosses(iCH4_) + Reaction

              IonSources(iCH3P_) = IonSources(iCH3P_) + Reaction
              NeutralSources(iH_) = NeutralSources(iH_) + Reaction

              ! ----------------------------------------------------------
              ! H2 + hv ==> Not H2
              ! ----------------------------------------------------------
              rr=EuvDissRateS(iLon,iLat,iAlt,iH2_,iBlock)

              Reaction = rr * Neutrals(iH2_)

              NeutralLosses(iH2_) = NeutralLosses(iH2_) + Reaction
!\
!-----------End Photochemistry-----------------------------------------------------+
!/
!
!
!
!
!\
!-----------Begin Bi-Molecular Neutral Chemistry-------------------------------------------+
!/
              ! -----------------------------------------------------------
              ! N2 + 1CH2 ==> 3CH2 + N2 
              ! -----------------------------------------------------------
              rr = (2.36e-14)*(tn)
              rr= rr*1.e-06  !cm3s-1-->m3s-1
              Reaction = rr * Neutrals(iN2_) * Neutrals(i1CH2_) 

              NeutralLosses(i1CH2_) = NeutralLosses(i1CH2_) + Reaction
              ! No Nitrogen Lost in this Reaction

              NeutralSources(i3CH2_) = NeutralSources(i3CH2_) + Reaction
!
              ! -----------------------------------------------------------
              ! CH4 + 1CH2 ==> 2*CH3 
              ! -----------------------------------------------------------
              rr = (6.0e-11)
              rr= rr*1.e-06  !cm3s-1-->m3s-1
              Reaction = rr * Neutrals(iCH4_) * Neutrals(i1CH2_) 

              NeutralLosses(i1CH2_) = NeutralLosses(i1CH2_) + Reaction
              NeutralLosses(iCH4_) = NeutralLosses(iCH4_) + Reaction

              NeutralSources(iCH3_) = NeutralSources(iCH3_) + 2.0*Reaction

              ! -----------------------------------------------------------
              ! 3CH2 + H ==> CH + H2
              ! -----------------------------------------------------------
              ! rt3CH2_H = exp(-370.0/tn)
              ! rr = (4.7e-10)*exp(-370.0/tn)

              rr = (4.7e-10)*rt3CH2_H

              rr= rr*1.e-06  !cm3s-1-->m3s-1
              Reaction = rr * Neutrals(i3CH2_) * Neutrals(iH_) 

              NeutralLosses(i3CH2_) = NeutralLosses(i3CH2_) + Reaction
              NeutralLosses(iH_) = NeutralLosses(iH_) + Reaction

              NeutralSources(iCH_) = NeutralSources(iCH_) + Reaction
              NeutralSources(iH2_) = NeutralSources(iH2_) + Reaction
!
              ! -----------------------------------------------------------
              ! CH + CH4 ==> C2H4 + H
              ! -----------------------------------------------------------
              ! rtCH_CH4 = (tn**(-1.04))*exp(-36.1/tn)
              ! rr = (3.96e-08)*(tn**(-1.04))*exp(-36.1/tn)

              rr = (3.96e-08)*rtCH_CH4

              rr= rr*1.e-06  !cm3s-1-->m3s-1
              Reaction = rr * Neutrals(iCH_) * Neutrals(iCH4_) 

              NeutralLosses(iCH_) = NeutralLosses(iCH_) + Reaction
              NeutralLosses(iCH4_) = NeutralLosses(iCH4_) + Reaction

              NeutralSources(iC2H4_) = NeutralSources(iC2H4_) + Reaction
              NeutralSources(iH_) = NeutralSources(iH_) + Reaction
  
              ! -----------------------------------------------------------
              ! 3CH2 + CH3 ==> C2H4 + H
              ! -----------------------------------------------------------
              rr = (7.00e-11)
              rr= rr*1.e-06  !cm3s-1-->m3s-1
              Reaction = rr * Neutrals(i3CH2_) * Neutrals(iCH3_) 

              NeutralLosses(i3CH2_) = NeutralLosses(i3CH2_) + Reaction
              NeutralLosses(iCH3_) = NeutralLosses(iCH3_) + Reaction

              NeutralSources(iC2H4_) = NeutralSources(iC2H4_) + Reaction
              NeutralSources(iH_) = NeutralSources(iH_) + Reaction
!
              ! -----------------------------------------------------------
              ! N(4S) + CH3 ==> H2CN + H
              ! -----------------------------------------------------------
              rr = (5.76e-11)

              rr= rr*1.e-06  !cm3s-1-->m3s-1
              Reaction = rr * Neutrals(iN4S_) * Neutrals(iCH3_) 

              NeutralLosses(iN4S_) = NeutralLosses(iN4S_) + Reaction
              NeutralLosses(iCH3_) = NeutralLosses(iCH3_) + Reaction

              NeutralSources(iH2CN_) = NeutralSources(iH2CN_) + Reaction
              NeutralSources(iH_) = NeutralSources(iH_) + Reaction
!
!              ! -----------------------------------------------------------
!              ! N(4S) + CH3 ==> HCN + H2
!              ! -----------------------------------------------------------
              rr = (6.00e-12)
              rr= rr*1.e-06  !cm3s-1-->m3s-1

              Reaction = rr * Neutrals(iN4S_) * Neutrals(iCH3_) 


              NeutralLosses(iN4S_) = NeutralLosses(iN4S_) + Reaction
              NeutralLosses(iCH3_) = NeutralLosses(iCH3_) + Reaction

              NeutralSources(iHCN_) = NeutralSources(iHCN_) + Reaction
              NeutralSources(iH2_) = NeutralSources(iH2_) + Reaction

           !   HCNChemicalSources(iLon,iLat,iAlt,1,iBlock) = Reaction
!
!             ! -----------------------------------------------------------
!             ! N(4S) + 3CH2 ==> HCN + H
!             ! -----------------------------------------------------------
              rr = (6.40e-12)
              rr= rr*1.e-06  !cm3s-1-->m3s-1
              Reaction = rr * Neutrals(iN4S_) * Neutrals(i3CH2_) 


             NeutralLosses(iN4S_) = NeutralLosses(iN4S_) + Reaction
             NeutralLosses(i3CH2_) = NeutralLosses(i3CH2_) + Reaction

             NeutralSources(iHCN_) = NeutralSources(iHCN_) + Reaction
             NeutralSources(iH_) = NeutralSources(iH_) + 2*Reaction
!
              ! -----------------------------------------------------------
              ! H2CN + H ==> HCN + H2
              ! -----------------------------------------------------------
              rr = (7.00e-11)
              rr= rr*1.e-06  !cm3s-1-->m3s-1
              Reaction = rr * Neutrals(iH2CN_) * Neutrals(iH_) 

              NeutralLosses(iH2CN_) = NeutralLosses(iH2CN_) + Reaction
              NeutralLosses(iH_) = NeutralLosses(iH_) + Reaction

              NeutralSources(iHCN_) = NeutralSources(iHCN_) + Reaction
              NeutralSources(iH2_) = NeutralSources(iH2_) + Reaction

            !  HCNChemicalSources(iLon,iLat,iAlt,2,iBlock) = Reaction
!
!              ! -----------------------------------------------------------
!              ! H2CN + N(4S) ==> HCN + NH
!              ! -----------------------------------------------------------
              rr = (3.98e-11)
              rr= rr*1.e-06  !cm3s-1-->m3s-1
              Reaction = rr * Neutrals(iH2CN_) * Neutrals(iN4S_) 


              NeutralLosses(iH2CN_) = NeutralLosses(iH2CN_) + Reaction
              NeutralLosses(iN4S_) = NeutralLosses(iN4S_) + Reaction
!
              NeutralSources(iHCN_) = NeutralSources(iHCN_) + Reaction
!
!\
!-----------End Bi-Molecular Neutral Chemistry-------------------------------------------+
!/
!\
!-----------Begin Bi-Molecular Ion-Neutral Chemistry-------------------------------------------+
!/
!
              ! -----------------------------------------------------------
              ! C2H5+ + HCN ==> C2H4 + H2CN+
              ! -----------------------------------------------------------
              ! +/- 20% error range
              ! -------------------]
              rr = (2.70e-09)*(0.80)
              rr= rr*1.e-06  !cm3s-1-->m3s-1
              Reaction = rr * Neutrals(iHCN_) * Ions(iC2H5P_) 

              NeutralLosses(iHCN_) = NeutralLosses(iHCN_) + Reaction
              IonLosses(iC2H5P_) = IonLosses(iC2H5P_) + Reaction

              NeutralSources(iC2H4_) = NeutralSources(iC2H4_) + Reaction
              IonSources(iH2CNP_) = IonSources(iH2CNP_) + Reaction

            !  HCNChemicalLosses(iLon,iLat,iAlt,1,iBlock) = Reaction
!
              ! -----------------------------------------------------------
              ! N+ + CH4 ==> CH3+ + NH 
              ! -----------------------------------------------------------
              ! +/- 15% error in reaction
              ! -------------------------]
              rr = (0.575e-09)
              rr= rr*1.e-06  !cm3s-1-->m3s-1
              Reaction = rr * Neutrals(iCH4_) * Ions(iNP_) 

              NeutralLosses(iCH4_) = NeutralLosses(iCH4_) + Reaction
              IonLosses(iNP_) = IonLosses(iNP_) + Reaction

              IonSources(iCH3P_) = IonSources(iCH3P_) + Reaction

              ! -----------------------------------------------------------
              ! N+ + CH4 ==> H2CN+ + H2 
              ! -----------------------------------------------------------
              rr = (0.414e-09)

              rr= rr*1.e-06  !cm3s-1-->m3s-1
              Reaction = rr * Neutrals(iCH4_) * Ions(iNP_) 

              NeutralLosses(iCH4_) = NeutralLosses(iCH4_) + Reaction
              IonLosses(iNP_) = IonLosses(iNP_) + Reaction
!
              NeutralSources(iH2_) = NeutralSources(iH2_) + Reaction
              IonSources(iH2CNP_) = IonSources(iH2CNP_) + Reaction
!
              ! -----------------------------------------------------------
              ! N2+ + CH4 ==> CH3+ + N2 + H 
              ! -----------------------------------------------------------
              ! +/- 15% error in this reaction
              ! -------------------------------]
              rr = (0.9804e-09)
              rr= rr*1.e-06  !cm3s-1-->m3s-1
              Reaction = rr * Neutrals(iCH4_) * Ions(iN2P_) 

              NeutralLosses(iCH4_) = NeutralLosses(iCH4_) + Reaction
              IonLosses(iN2P_) = IonLosses(iN2P_) + Reaction

              NeutralSources(iN2_) = NeutralSources(iN2_) + Reaction
              NeutralSources(iH_) = NeutralSources(iH_) + Reaction
              IonSources(iCH3P_) = IonSources(iCH3P_) + Reaction

              ! Isotope Sources
              NeutralSources(i15N2_) = NeutralSources(i15N2_) + &
                      0.5*(Neutrals(i15N2_)/Neutrals(iN2_))*Reaction
!
              ! -----------------------------------------------------------
              ! CH3+ + CH4 ==> C2H5+ + H2 
              ! -----------------------------------------------------------
              ! +/- 20% error in this reaction
              ! ------------------------------]
              rr = (1.10e-09)
              rr= rr*1.e-06  !cm3s-1-->m3s-1
              Reaction = rr * Neutrals(iCH4_) * Ions(iCH3P_) 

              NeutralLosses(iCH4_) = NeutralLosses(iCH4_) + Reaction
              IonLosses(iCH3P_) = IonLosses(iCH3P_) + Reaction

              NeutralSources(iH2_) = NeutralSources(iH2_) + Reaction
              IonSources(iC2H5P_) = IonSources(iC2H5P_) + Reaction
!
!\
!-----------END Bi-Molecular Ion-Neutral Chemistry-------------------------------------------+
!/
!\
!-----------Begin Electron Recombination Chemistry-------------------------------------------+
!/
              ! -----------------------------------------------------------
              ! H2CN+ + e- ==> HCN + H 
              ! -----------------------------------------------------------
              ! rt300te = sqrt(300.0/te)
              ! rr = (6.40e-07)*sqrt(300.0/te)

              rr = (6.40e-07)*rt300te
              rr= rr*1.e-06  !cm3s-1-->m3s-1

              Reaction = rr * Ions(ie_) * Ions(iH2CNP_) 

              IonLosses(iH2CNP_) = IonLosses(iH2CNP_) + Reaction

              NeutralSources(iH_) = NeutralSources(iH_) + Reaction
              NeutralSources(iHCN_) = NeutralSources(iHCN_) + Reaction
!
            !  HCNChemicalSources(iLon,iLat,iAlt,3,iBlock) = Reaction
              ! -----------------------------------------------------------
              ! C2H5+ + e- ==> C2H4 + H 
              ! -----------------------------------------------------------
              ! rt300te = sqrt(300.0/te)
              ! rr = (0.072e-06)*sqrt(300.0/te)

              rr = (0.072e-06)*rt300te

              rr= rr*1.e-06  !cm3s-1-->m3s-1
              Reaction = rr * Ions(ie_) * Ions(iC2H5P_) 

              IonLosses(iC2H5P_) = IonLosses(iC2H5P_) + Reaction
!
              NeutralSources(iH_) = NeutralSources(iH_) + Reaction
              NeutralSources(iC2H4_) = NeutralSources(iC2H4_) + Reaction

!
!\
!-----------End Electron Recombination Chemistry-------------------------------------------+
!/
!
! \ 
! Begin Special C2H4 Losses 
!              ! -----------------------------------------------------------
!              ! C2H4 + CH ==> (stuff) + H 
!              !    And 
!              ! C2H4 + CH ==> (other stuff) + H 
!              ! -----------------------------------------------------------
!              ! rtC2H4_CH = (tn**(-.546))*(exp(-29.6/tn))
!              ! rr = (3.87e-09)*(tn**(-.546))*(exp(-29.6/tn))
!
!              rr = (3.87e-09)*rtC2H4_CH
!
!              rr= rr*1.e-06  !cm3s-1-->m3s-1
!
!              Reaction = rr * Neutrals(iC2H4_) * Neutrals(iCH_) 
!
!              NeutralLosses(iC2H4_) = NeutralLosses(iC2H4_) + 2.0*Reaction
!              NeutralLosses(iCH_) = NeutralLosses(iCH_) + 2.0*Reaction
!!
!              NeutralSources(iH_) = NeutralSources(iH_) + 2.0*Reaction
!
!              StuffSources = StuffSources + 2.0*Reaction*(Mass(iC2H4_) + 12.0*AMU)
!
!              ! -----------------------------------------------------------
!              ! CH3+ + C2H4 ==> (stuff) + H2 
!              ! -----------------------------------------------------------
!              rr = (.5406e-09)
!              rr= rr*1.e-06  !cm3s-1-->m3s-1
!
!              Reaction = rr * Neutrals(iC2H4_) * Ions(iCH3P_) 
!
!              NeutralLosses(iC2H4_) = NeutralLosses(iC2H4_) + Reaction
!              IonLosses(iCH3P_) = IonLosses(iCH3P_) + Reaction
!
!              NeutralSources(iH2_) = NeutralSources(iH2_) + Reaction
!
!              StuffSources = StuffSources + Reaction*(Mass(iC2H4_) + Mass(iCH_))
!!
!! /
!              ! -----------------------------------------------------------
!              ! CH3+ + C2H4 ==> (stuff) + CH4 
!              ! -----------------------------------------------------------
!              rr = (.4879e-09)
!              rr= rr*1.e-06  !cm3s-1-->m3s-1
!
!              Reaction = rr * Neutrals(iC2H4_) * Ions(iCH3P_) 
!
!              NeutralLosses(iC2H4_) = NeutralLosses(iC2H4_) + Reaction
!              IonLosses(iCH3P_) = IonLosses(iCH3P_) + Reaction
!
!              NeutralSources(iCH4_) = NeutralSources(iCH4_) + Reaction
!
!              StuffSources = StuffSources + Reaction*(Mass(iC2H4_) - Mass(iH_))
!
!              ! -----------------------------------------------------------
!              ! C2H5+ + C2H4 ==> (stuff) + CH4 
!              ! -----------------------------------------------------------
!              rr = (.355e-09)
!              rr= rr*1.e-06  !cm3s-1-->m3s-1
!
!              Reaction = rr * Neutrals(iC2H4_) * Ions(iC2H5P_) 
!
!              NeutralLosses(iC2H4_) = NeutralLosses(iC2H4_) + Reaction
!              IonLosses(iC2H5P_) = IonLosses(iC2H5P_) + Reaction
!
!              NeutralSources(iCH4_) = NeutralSources(iCH4_) + Reaction
!
!              StuffSources = StuffSources + Reaction*(Mass(iC2H5P_) + 12.0*AMU)
!
!              ! -----------------------------------------------------------
!              ! N+ + C2H4 ==> H2CN+ + 3CH2 
!              ! -----------------------------------------------------------
!              rr = (.195e-09)
!              rr= rr*1.e-06  !cm3s-1-->m3s-1
!
!              Reaction = rr * Neutrals(iC2H4_) * Ions(iNP_) 
!
!              NeutralLosses(iC2H4_) = NeutralLosses(iC2H4_) + Reaction
!              IonLosses(iNP_) = IonLosses(iNP_) + Reaction
!
!              NeutralSources(i3CH2_) = NeutralSources(i3CH2_) + Reaction
!              IonSources(iH2CNP_) = IonSources(iH2CNP_) + Reaction
!!
!              ! -----------------------------------------------------------
!              ! N+ + C2H4 ==> Stuff + N(4S) 
!              ! -----------------------------------------------------------
!              rr = (.455e-09)
!              rr= rr*1.e-06  !cm3s-1-->m3s-1
!
!              Reaction = rr * Neutrals(iC2H4_) * Ions(iNP_) 
!
!              NeutralLosses(iC2H4_) = NeutralLosses(iC2H4_) + Reaction
!              IonLosses(iNP_) = IonLosses(iNP_) + Reaction
!
!              NeutralSources(iN4S_) = NeutralSources(iN4S_) + Reaction
!
!              StuffSources = StuffSources + Reaction*(Mass(iC2H4_))
!
!              ! -----------------------------------------------------------
!              ! N2+ + C2H4 ==> Stuff + HCN + H2 
!              ! -----------------------------------------------------------
!              rr = (.130e-09)*1.15
!!              rr= rr*1.e-06  !cm3s-1-->m3s-1
!
!              Reaction = rr * Neutrals(iC2H4_) * Ions(iN2P_) 
!
!              NeutralLosses(iC2H4_) = NeutralLosses(iC2H4_) + Reaction
!              IonLosses(iN2P_) = IonLosses(iN2P_) + Reaction
!
!              NeutralSources(iHCN_) = NeutralSources(iHCN_) + Reaction
!              NeutralSources(iH2_) = NeutralSources(iH2_) + Reaction
!
!              StuffSources = StuffSources + Reaction*(AMU + 12.0*AMU + 14.0*AMU  )
!
!           !   HCNChemicalSources(iLon,iLat,iAlt,4,iBlock) = Reaction
!              ! -----------------------------------------------------------
!              ! N2+ + C2H4 ==> H2CN+ + HCN + H 
!              ! -----------------------------------------------------------
!              rr = (.130e-09)*1.15
!              rr= rr*1.e-06  !cm3s-1-->m3s-1
!
!              Reaction = rr * Neutrals(iC2H4_) * Ions(iN2P_) 
!
!              NeutralLosses(iC2H4_) = NeutralLosses(iC2H4_) + Reaction
!              IonLosses(iN2P_) = IonLosses(iN2P_) + Reaction
!
!              IonSources(iH2CNP_) = IonSources(iH2CNP_) + Reaction
!              NeutralSources(iHCN_) = NeutralSources(iHCN_) + Reaction
!              NeutralSources(iH_) = NeutralSources(iH_) + Reaction
!
!           !   HCNChemicalSources(iLon,iLat,iAlt,5,iBlock) = Reaction
!
!! ----End Special C2H4 Loss Section ---------------------------->>>
!!/
!

!
! \
! CH4 Aerosol Loss Rate
!              ! ----------------------------------------------------------
!              ! Ar + Aerosol ==> Lost CH4
!              ! ----------------------------------------------------------
!
!                Reaction = ArAerosolLosses 
!
!                NeutralLosses(iAr_) = NeutralLosses(iAr_) + Reaction
!!              ! ----------------------------------------------------------
!              ! CH4 + Aerosol ==> Lost CH4
!              ! ----------------------------------------------------------
!
!                Reaction = CH4AerosolLosses 
!
!                NeutralLosses(i13CH4_) = &
!                      (0.9*NeutralLosses(iCH4_) + 0.9*Reaction)*(Neutrals(i13CH4_)/Neutrals(iCH4_)) 
!
!                NeutralLosses(iCH4_) = NeutralLosses(iCH4_) + Reaction
!
!!              ! ----------------------------------------------------------
!!              ! H2 + Aerosol ==> Lost H2
!!              ! ----------------------------------------------------------
!!
!                Reaction = H2AerosolLosses 
!
!                NeutralLosses(iH2_) = &
!                      NeutralLosses(iH2_) + Reaction
!
!!              ! ----------------------------------------------------------
!!              ! N2 + Aerosol ==> Lost N2
!!              ! ----------------------------------------------------------
!!
!                Reaction = N2AerosolLosses 
!
!                NeutralLosses(iN2_) = &
!                      NeutralLosses(iN2_) + Reaction
!
!!              ! ----------------------------------------------------------
!!              ! H + Aerosol ==> Produced H2
!!              ! ----------------------------------------------------------
!!
!                Reaction = H2AerosolSources 
!
!                NeutralSources(iH2_) = &
!                      NeutralSources(iH2_) + Reaction
!
!                NeutralLosses(iH_) = &
!                      NeutralLosses(iH_) + 2.0*Reaction
!
!\
!\
!------------------------------------------------------+
!-----------End Total Chemistry-------------------------------------------+
!------------------------------------------------------+
!/


!\
! Special Section:----------------
! Chemical Sources for 15N2_, 13CH4_

!                NeutralLosses(i13CH4_) = &
!                      (0.9*NeutralLosses(iCH4_))*(Neutrals(i13CH4_)/Neutrals(iCH4_)) 

             NeutralSources(i13CH4_) = 0.9*NeutralSources(iCH4_)*(Neutrals(i13CH4_)/Neutrals(iCH4_))


!             NeutralSources(i15N2_) = 0.0
!             NeutralLosses(i15N2_) = 0.0





           !   HCNChemicalSources(iLon,iLat,iAlt,6,iBlock) = NeutralSources(iHCN_)
           !   HCNChemicalLosses(iLon,iLat,iAlt,2,iBlock) = NeutralLosses(iHCN_)

!\
! Freezing HCN Chemistry for now (4-6-2006)
!
!             NeutralSources(iHCN_) = 0.0
!             NeutralLosses(iHCN_) = 0.0
!
!             NeutralSources(iN2_) = 0.0
!             NeutralLosses(iN2_) = 0.0
!
!             NeutralSources(iCH4_) = 0.0
!             NeutralLosses(iCH4_) = 0.0
!
             NeutralSources(iH2_) = 0.0
             NeutralLosses(iH2_) = 0.0
!\
!  Temporarily Freezing things!!
!/
!
!              NeutralSources(1:nSpeciesTotal) = 0.0
!              NeutralLosses(1:nSpeciesTotal) = 0.0
!
!  Freezing Ion Chemistry from above !!
!  Full Ion Equilibrium Chemistry Done below 
!              NeutralSources(i3CH2_) = 0.0
!              NeutralLosses(i3CH2_) = 0.0
!
!              NeutralSources(iHCN_) = 0.0
!              NeutralLosses(iHCN_) = 0.0

!              NeutralSources(iH2_) = 0.0
!              NeutralLosses(iH2_) = 0.0
!
!              IonSources(1:nIons) = 0.0
!              IonLosses(1:nIons) = 0.0
!
!              IonSources(1) = 0.0
!              IonLosses(1) = 0.0
!              IonSources(3:nIons) = 0.0
!              IonLosses(3:nIons) = 0.0
!
!/
              
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

              tli = DtSub * IonLosses
              tsi = DtSub * IonSources + Ions

              do iIon = 1, nIons-1
                 do while (tsi(iIon)-tli(iIon) < 0.0 .and. DtSub > 1.0e-2)
                    if (tsi(iIon)-tli(iIon) < 0.0 .and. Ions(iIon) < 1.0e7) then
                       IonLosses(iIon) = &
                            (IonSources(iIon) + Ions(iIon)/DtSub)*0.9
                    else
                       DtSub = DtSub/2.0
                    endif
                    tli(iIon) = DtSub * IonLosses(iIon)
                    tsi(iIon) = DtSub * IonSources(iIon) + Ions(iIon)
                 enddo
              enddo

              !---- Neutrals

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
              
              tln = DtSub * NeutralLosses
              tsn = DtSub * NeutralSources + 0.1*Neutrals
              do while (minval(tsn-tln) < 0.0)
                 DtSub = DtSub/2.0
                 tln = DtSub * NeutralLosses
                 tsn = DtSub * NeutralSources + 0.1*Neutrals
              enddo

              Ions(nIons) = 0.0
!              !if (ialt .eq. 10) write(*,*) "NO sources and losses: ",NeutralSources(iNO_), &
!                   NeutralLosses(iNO_)
              do iIon = 1, nIons-1

                 if (Ions(iIon) + &
                      (IonSources(iIon) - IonLosses(iIon)) * DtSub < 0.0) then
                    !!!!!! Solve Steady-State !!!!!!!
                    Ions(iIon) = IonSources(iIon)*Ions(iIon)/IonLosses(iIon)
                 else
                    Ions(iIon) = Ions(iIon) + &
                         (IonSources(iIon) - IonLosses(iIon)) * DtSub
                 endif

!                 Ions(iIon) = max(0.01,Ions(iIon))
                 
                 ! sum for e-
                 Ions(nIons) = Ions(nIons) + Ions(iIon)

                 if (Ions(iIon) < 0.0) then
                    write(*,*) "Negative Ion Density : ", &
                         iIon, iLon, iLat, iAlt, &
                         Ions(iIon), &
                         IonSources(iIon), IonLosses(iIon)
                 endif
              enddo

              do iNeutral = 1, nSpeciesTotal
                 Neutrals(iNeutral) = &
                      Neutrals(iNeutral) + &
                      (NeutralSources(iNeutral) - NeutralLosses(iNeutral)) * &
                      DtSub
                 
                 NeutralSourcesTotal(iNeutral,iAlt) = NeutralSourcesTotal(iNeutral,iAlt) + &
                      NeutralSources(iNeutral) * DtSub

                 NeutralLossesTotal(iNeutral,iAlt) = NeutralLossesTotal(iNeutral,iAlt) + &
                      NeutralLosses(iNeutral) * DtSub
                 

                 if (Neutrals(iNeutral) < 0.0) then
                    write(*,*) "Negative Neutral Density : ", &
                         iNeutral, iLon, iLat, iAlt, DtSub, &
                         Neutrals(iNeutral), &
                         NeutralSources(iNeutral), NeutralLosses(iNeutral)
                 endif


              enddo

              ChemicalHeatingRate(iLon,iLat,iAlt) = &
                   ChemicalHeatingRate(iLon,iLat,iAlt) + &
                   ChemicalHeatingSub * DtSub

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
                    write(*,*) "Neutral Source/Loss : ", &
                         iNeutral, NeutralSources(iNeutral), &
                         NeutralLosses(iNeutral)
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
