!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine calc_chemistry(iBlock)

  use ModSizeGitm
  use ModGITM
  use ModPlanet
  use ModRates
  use ModEUV
  use ModSources
  use ModInputs, only: iDebugLevel, UseIonChemistry, UseNeutralChemistry
  use ModConstants

  implicit none

  integer, intent(in) :: iBlock

  real :: IonSources(nIons), NeutralSources(nSpeciesTotal)
  real :: IonLosses(nIons), NeutralLosses(nSpeciesTotal)
  real :: DtSub, DtOld, DtTotal, DtMin, DtAve, Source, Reaction
  real :: Ions(nIons), Neutrals(nSpeciesTotal), tsi(nIons), tli(nIons)
  real :: tln(nSpeciesTotal), tsn(nSpeciesTotal)

  integer :: iLon, iLat, iAlt, iIon, nIters, iNeutral
  
  real :: ChemicalHeatingSub, rr
  real :: Emission(nEmissions), EmissionTotal(nEmissions)

  logical :: UseNeutralConstituent(nSpeciesTotal)
  logical :: UseIonConstituent(nIons)
  !---------------------------------------------------------------------------

  UseNeutralConstituent = .true.
  UseIonConstituent     = .true.

  DtMin = Dt

  call report("Chemistry",2)
  call start_timing("calc_chemistry")

  DtAve = 0.0

  nIters=0

  do iLon = 1, nLons
     do iLat = 1, nLats
        do iAlt = 1, nAlts

           NeutralSourcesTotal = 0.0
           NeutralLossesTotal = 0.0

           DtTotal = 0.0
           EmissionTotal = 0.0

           Ions = IDensityS(iLon,iLat,iAlt,:,iBlock)

           Neutrals = NDensityS(iLon,iLat,iAlt,:,iBlock)

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

              rr=EuvDissRateS(iLon,iLat,iAlt,iH2_,iBlock)
  
              Reaction = rr * &
                   Neutrals(iH2_)
  
              NeutralLosses(iH2_) = NeutralLosses(iH2_) + Reaction
              NeutralSources(iH_) = NeutralSources(iH_) + 2*Reaction

              ! Solar EUV
  
              Reaction = EuvIonRateS(iLon,iLat,iAlt,iH2P_,iBlock) * &
                   Neutrals(iH2_)
  
              IonSources(iH2P_)   = IonSources(iH2P_)   + Reaction
              NeutralLosses(iH2_) = NeutralLosses(iH2_) + Reaction
              
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
                 
                 NeutralSourcesTotal(ialt,iNeutral) = &
                      NeutralSourcesTotal(ialt,iNeutral) + &
                      NeutralSources(iNeutral) * DtSub

                 NeutralLossesTotal(ialt,iNeutral) = &
                      NeutralLossesTotal(ialt,iNeutral) + &
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
 
  call end_timing("calc_chemistry")

end subroutine calc_chemistry
