!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine calc_chemistry(iBlock)

  use ModSizeGitm
  use ModGITM
  use ModPlanet
  use ModRates
  use ModEUV
  use ModSources
  use ModInputs, only: iDebugLevel, UseIonChemistry, UseNeutralChemistry,f107
  use ModConstants

  implicit none

  integer, intent(in) :: iBlock

  real :: IonSources(nIons), NeutralSources(nSpeciesTotal)
  real :: IonLosses(nIons), NeutralLosses(nSpeciesTotal)
  real :: DtSub, DtOld, DtTotal, DtMin, DtAve, Source, Reaction
  real :: tn, dtsubtmp, losstmp, dentmp
  real :: Ions(nIons), Neutrals(nSpeciesTotal)

  integer :: iLon, iLat, iAlt, iIon, nIters, iDtReducer, iNeutral

  real :: ChemicalHeatingSub
  real :: Emission(nEmissions), EmissionTotal(nEmissions)
  real :: te300, rr01, rr02, rr03, rr11, rr12, nn

  logical :: UseNeutralConstituent(nSpeciesTotal)
  logical :: UseIonConstituent(nIons)

  UseNeutralConstituent = .true.
  UseIonConstituent     = .true.

  DtMin = Dt

  if (.not.UseIonChemistry) return

  call report("Chemistry",2)
  call start_timing("calc_chemistry")

  DtAve = 0.0

  nIters=0

  do iLon = 1, nLons
     do iLat = 1, nLats
        do iAlt = 1, nAlts

           tn = Temperature(iLon,iLat,iAlt,iBlock)/TempUnit(iLon,iLat,iAlt)

           !! CHANGE
           te300 = 300.0 / eTemperature(iLon,iLat,iAlt,iBlock)

           rr01 = 1.5e-29 / (tn ^ (1.0/3.0))
           rr02 = 3.73e-20 * tn^3/exp(4406.0/tn)
           rr03 = 6.40e-25 / (tn^2 * exp(1200.0/tn))

           rr11 = 4.8e-8 * te300^0.5
           rr12 = 6.2e-8 * te300^0.5

           DtTotal = 0.0
           EmissionTotal = 0.0

           Ions = IDensityS(iLon,iLat,iAlt,:,iBlock)
           Neutrals = NDensityS(iLon,iLat,iAlt,:,iBlock)
           nn = sum(Neutrals)
 
           do while (DtTotal < Dt)

              ChemicalHeatingSub = 0.0
              Emission = 0.0

              DtSub = Dt - DtTotal

              IonSources = 0.0
              NeutralSources = 0.0
              IonLosses  = 0.0
              NeutralLosses = 0.0

              ! Solar EUV

              ! ----------------------------------------------------------
              ! H+
              ! ----------------------------------------------------------

              ! Solar EUV

              Reaction = EuvIonRateS(iLon,iLat,iAlt,iHP_,iBlock) * &
                   Neutrals(iH_)

              IonSources(iHP_)   = IonSources(iHP_)   + Reaction
              NeutralLosses(iH_) = NeutralLosses(iH_) + Reaction

              ! 2H + M -> H2 + M

              Reaction = &
                   rr01 * 2*Neutrals(iH_) * nn

              NeutralSources(iH2_) = NeutralSources(iH2_) + Reaction
              NeutralLosses(iH_)   = NeutralLosses(iH_)   + 2*Reaction

              ! H + CH4 -> CH3 + H2

              Reaction = &
                   rr02 * Neutrals(iH_) * Neutrals(iCH4_)

              NeutralSources(iH2_)  = NeutralSources(iH2_) + Reaction
              NeutralSources(iCH3_) = NeutralSources(iCH3_) + Reaction
              NeutralLosses(iH_)    = NeutralLosses(iH_) + Reaction
              NeutralLosses(iCH4_)  = NeutralLosses(iCH4_) + Reaction

              ! H + C2H2 + M -> C2H3 + M

              Reaction = &
                   rr03 * Neutrals(iH_) * Neutrals(iCH4)

              NeutralSources(iC2H3_) = NeutralSources(iC2H3_) + Reaction
              NeutralLosses(iH_)     = NeutralLosses(iH_) + Reaction
              NeutralLosses(iC2H2_)  = NeutralLosses(iC2H2_) + Reaction

              ! H2+ + H -> H+ + H2

              Reaction = &
                   RrTempInd(iRrR4_) * Ions(iH2P_) * Neutrals(iH_)

              IonSources(iHP_)   = IonSources(iHP_)   + Reaction
              NeutralSources(iH2_) = NeutralSources(iH2_) + Reaction
              IonLosses(iH2P_)    = IonLosses(iH2P_)  + Reaction
              NeutralLosses(iH_)  = NeutralLosses(iH_) + Reaction

              ! H2+ + H2 -> H3+ + H

              Reaction = &
                   RrTempInd(iRrR5_) * Ions(iH2P_) * Neutrals(iH2_)

              IonSources(iH3P_)   = IonSources(iH3P_)   + Reaction
              NeutralSources(iH_) = NeutralSources(iH_) + Reaction
              IonLosses(iH2P_)    = IonLosses(iH2P_)    + Reaction
              NeutralLosses(iH2_) = NeutralLosses(iH2_) + Reaction

              ! H2+ + He -> HeH+ + H

              Reaction = &
                   RrTempInd(iRrR6_) * Ions(iH2P_) * Neutrals(iHe_)

              IonSources(iCH5P_)   = IonSources(iCH5P_)   + Reaction
              NeutralSources(iH_) = NeutralSources(iH_) + Reaction
              IonLosses(iH2P_)    = IonLosses(iH2P_)    + Reaction
              NeutralLosses(iH2_) = NeutralLosses(iH2_) + Reaction


!!!
!!!              !CO2+ + NO -> NO+ + CO2
!!!
!!!              Reaction = &
!!!                   RrTempInd(iRrR5_) * &
!!!                   Ions(iCO2P_) * &
!!!                   Neutrals(iNO_)
!!!
!!!              IonSources(iNOP_)     = IonSources(iNOP_)   + Reaction
!!!              NeutralSources(iCO2_) = NeutralSources(iCO2_) + Reaction
!!!              IonLosses(iCO2P_)     = IonLosses(iCO2P_)  + Reaction
!!!              NeutralLosses(iCO2_)  = NeutralLosses(iCO2_) + Reaction
!!!
!!!!!              !CO+ + NO -> NO+ + CO
!!!!!
!!!!!              Reaction = &
!!!!!                   RrTempInd(iRrR7_) * &
!!!!!                   Ions(iCOP_) * &
!!!!!                   Neutrals(iNO_)
!!!!!
!!!!!              IonSources(iNOP_)     = IonSources(iNOP_)   + Reaction
!!!!!              NeutralSources(iCO_) = NeutralSources(iCO_) + Reaction
!!!!!              IonLosses(iCOP_)     = IonLosses(iCOP_)  + Reaction
!!!!!              NeutralLosses(iNO_)  = NeutralLosses(iNO_) + Reaction
!!!
!!!!!              ChemicalHeatingSub = &
!!!!!                   ChemicalHeatingSub + Reaction * 1.33

              ! ----------------------------------------------------------
              ! CO2+
              ! ----------------------------------------------------------

              ! Solar EUV

              Reaction = EuvIonRateS(iLon,iLat,iAlt,iCO2P_,iBlock) * &
                   Neutrals(iCO2_)

              IonSources(iCO2P_)   = IonSources(iCO2P_)   + Reaction
              NeutralLosses(iCO2_) = NeutralLosses(iCO2_) + Reaction



!!              ! Aurora
!!
!!              Reaction = AuroralIonRateS(iLon,iLat,iAlt,iN2_, iBlock) + &
!!                   IonPrecipIonRateS(iLon,iLat,iAlt,iN2_, iBlock)
!!
!!              IonSources(iN2P_)   = IonSources(iN2P_) + Reaction
!!              NeutralLosses(iN2_) = NeutralLosses(iN2_) + Reaction
!!
              


              do iIon = 1, nIons-1

                 if (UseIonChemistry .and. &
                      UseIonConstituent(iIon)) then

                    if (IonLosses(iIon) > IonSources(iIon) .and. &
                         IonSources(iIon)>0.0) then
                       if (Ions(iIon) > 0.01) then
                          dtOld = DtSub
                          dtSub = &
                               min(0.25 * &
                               (IonSources(iIon) + &
                               Ions(iIon))/ &
                               (abs(IonLosses(iIon))+0.1), DtOld)
                          if (DtSub < DtOld) iDtReducer = iIon
                       else
                          IonSources(iIon) = 0.0
                          IonLosses(iIon) = 0.0
                       endif
                    endif
                 else
                    IonSources(iIon) = 0.0
                    IonLosses(iIon) = 0.0
                 endif
              enddo

              do iNeutral = 1, nSpeciesTotal

                 if (UseNeutralChemistry .and. &
                      UseNeutralConstituent(iNeutral)) then

                    if (NeutralLosses(iNeutral) > &
                         NeutralSources(iNeutral)) then
                       if (Neutrals(iNeutral)>0.01) then
                          dtOld = DtSub
                          dtSub = &
                               min(0.25 * &
                               Neutrals(iNeutral)/ &
                               (abs(NeutralSources(iNeutral) - &
                               NeutralLosses(iNeutral))+0.1), DtOld)
                          if (DtSub < DtOld) iDtReducer = nIons + iNeutral
                       else
                          NeutralSources(iNeutral) = 0.0
                          NeutralLosses(iNeutral) = 0.0
                       endif
                    endif
                 else
                    NeutralSources(iNeutral) = 0.0
                    NeutralLosses(iNeutral) = 0.0
                 endif
              enddo

              Ions(nIons) = 0.0
              do iIon = 1, nIons-1
                 Ions(iIon) = Ions(iIon) + &
                      (IonSources(iIon) - IonLosses(iIon)) * DtSub
                 Ions(iIon) = max(0.01,Ions(iIon))
                 
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
              enddo

              ChemicalHeatingRate(iLon,iLat,iAlt) = &
                   ChemicalHeatingRate(iLon,iLat,iAlt) + &
                   ChemicalHeatingSub * DtSub

              EmissionTotal = EmissionTotal + Emission(:)*DtSub

              DtTotal = DtTotal + DtSub

              if (DtSub < DtMin) DtMin = DtSub

              if (DtSub < 1.0e-9 .and. abs(DtTotal-Dt) > DtSub) then
                 write(*,*) "Chemistry is too fast!!", DtSub
                 if (iDtReducer < nIons) then
                    write(*,*) "Ion Constituent : ", &
                         iDtReducer, iLon, iLat, iAlt,&
                         Ions(iDtReducer), &
                         IonSources(iDtReducer), IonLosses(iDtReducer)
                 else
                    iDtReducer = iDtReducer - nIons
                    write(*,*) "Neutral Constituent : ", &
                         iDtReducer, iLon, iLat, iAlt,&
                         Ions(iDtReducer), &
                         IonSources(iDtReducer), IonLosses(iDtReducer)
                 endif

                 call stop_gitm("Ion Chemistry is too fast!!")
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
        write(*,*) "====> Max Ion Density: ", iIon, &
             maxval(IDensityS(1:nLons,1:nLats,(nAlts*4)/5,iIon,iBlock))
     enddo
  endif

  if (iDebugLevel > 2) &
       write(*,*) "===> Average Dt for this timestep : ", &
       (Dt*nLats*nLons*nAlts)/nIters

  call end_timing("calc_chemistry")

end subroutine calc_chemistry
