subroutine calc_chemistry(iBlock)

!  No override of CO2 densities; explicit calculation

  use ModGITM
  use ModInputs, only: iDebugLevel, UseIonChemistry, UseNeutralChemistry,f107,f107a, &
       useImplicitChemistry
  use ModConstants
  use ModSources
  use ModChemistry
  use ModGITMImplicit
  use ModTime, only : iStep
  use ModPlanet, only : ialtminiono
  use ModUserGITM

  implicit none

  integer, intent(in) :: iBlock

  real  :: IonSources(nIons),IonLosses(nIons),nSources(nSpeciesTotal),ISources(nIons)
  real  :: NeutralSources(nSpeciesTotal), NeutralLosses(nSpeciesTotal)
  real  :: ChemicalHeatingSub,Emission(nEmissions)

  integer :: iLon,iLat,iAlt,niters,iIon, iNeutral, ivar,nImplicitSourceCalls,iminiono

  real :: dttotal, dtsub, dtMin
  real :: tli(nIons), tsi(nIons), tln(nSpeciesTotal), tsn(nSpeciesTotal)

  real :: EmissionTotals(nEmissions), F

  logical :: doImplicit = .true.
  !---------------------------------------------------------------------------
  UseNeutralConstituent = .true.
  UseIonConstituent     = .true.

  nImplicitSourceCalls = nSpeciesAll * 2 + 1

  call report("Chemistry",2)
  call start_timing("calc_chemistry")
  if (iDebugLevel > 3) then
     do iIon = 1, nIons
        write(*,*) "====> start calc_chemistry: Max Ion Density: ", iIon, &
             maxval(IDensityS(1:nLons,1:nLats,(nAlts*4)/5,iIon,iBlock))
     enddo
  endif

  !   if (istep .lt. 10000) then
  !      useimplicitchemistry = .false.
  !   else 
  !      useimplicitchemistry = .true.
  !   endif
!  neutralsourcestotal = 0.0
!  neutrallossestotal = 0.0

  do iLon = 1, nLons
     do iLat = 1, nLats

        iMinIono = iAltMinIono(iLon,iLat,iBlock)

        do ialt = iminiono, nAlts

           dtTotal = 0.0
           DtMin = Dt
           DtSub = Dt - DtTotal

           Ions = IDensityS(iLon,iLat,iAlt,1:nIons,iBlock)
           Neutrals = NDensityS(iLon,iLat,iAlt,:,iBlock)

           call calc_reaction_rates(iLon,iLat,iAlt,iBlock)
           call calc_chemical_sources(iLon,iLat,iAlt,iBlock,IonSources,IonLosses,NeutralSources, &
                NeutralLosses,ChemicalHeatingSub,Emission)

           !        write(*,*) dtsub
           call calc_dtsub(IonSources,IonLosses,NeutralSources,NeutralLosses,dtSub)

           !          write(*,*)ialt, "predicted: ",dt/dtsub,dtsub,dt,nimplicitsourcecalls
           if (.05*Dt/DtSub > nImplicitSourceCalls .and. useImplicitChemistry) then
              doImplicit = .true.
           else
              doImplicit = .false.
           endif


           if (doImplicit) then
              !write(*,*)ialt," implicit"
              where (SpeciesDensity(:,:,:,:,iBlock) .lt. 0) &
                   SpeciesDensity(:,:,:,:,iBlock) = 0.0

              where (SpeciesDensityOld(:,:,:,:,iBlock) .lt. 0) &
                   SpeciesDensityOld(:,:,:,:,iBlock) = 0.0


              EmissionTotals = 0.0

              call update_point_implicit(iLon,iLat,iAlt,iBlock,nSources,ISources,ChemicalHeatingSub)

              do iVar = 0, nSpeciesTotal
                 if (.not. UseNeutralConstituent(ivar)) nSources(ivar) = 0
              enddo

              do iVar = 0, nIons
                 if (.not. UseIonConstituent(ivar)) iSources(ivar) = 0
              enddo

              NDensityS(iLon,iLat,iAlt,:,iBlock) =  NDensityS(iLon,iLat,iAlt,:,iBlock) + nSources
              IDensityS(iLon,iLat,iAlt,:,iBlock)   =  IDensityS(iLon,iLat,iAlt,:,iBlock) + iSources 

              IDensityS(iLon,iLat,iAlt,nIons,iBlock) = 0.0
              do iIon = 1, nIons - 1
                 IDensityS(iLon,iLat,iAlt,nIons,iBlock) = IDensityS(iLon,iLat,iAlt,nIons,iBlock) + &
                      IDensityS(iLon,iLat,iAlt,iIon,iBlock)
              enddo

              ! end if implicit --------------------------
           else 
              !write(*,*) ialt,'Explicit'

              niters = 0

              do iIon = 1, nIons-1
                 if (IonSources(iion) > 100.0*ions(iion) .and. IonLosses(iion) > 100.0 * ions(iion)) then
                    F = ions(iion)/IonSources(iion) * 10.0
                    IonLosses(iion) = f * IonLosses(iion)
                    IonSources(iion) = f * IonSources(iion)
                 endif
              enddo

              Ions(nIons) = 0.0
              do iIon = 1, nIons - 1
                 Ions(iIon) = Ions(iIon) + &
                      (IonSources(iIon) - IonLosses(iIon))*DtSub

                 ! sum for e-
                 Ions(nIons) = Ions(nIons) + Ions(iIon)

                 if (Ions(iIon) < 0.0) then
                    write(*,*) "Negative Ion Density : ", &
                         iIon, iLon, iLat, iAlt, &
                         Ions(iIon), &
                         IonSources(iIon), IonLosses(iIon)
                 endif
              enddo

              userdata1d(1,1,ialt,25) = 0
              userdata1d(1,1,ialt,25) = (NeutralSources(1)-NeutralLosses(1)) * DtSub
              do iNeutral = 1, nSpeciesTotal
                 Neutrals(iNeutral) = &
                      Neutrals(iNeutral) + &
                      (NeutralSources(iNeutral) - NeutralLosses(iNeutral)) * &
                      DtSub

!                 NeutralSourcesTotal(iLon,iLat,iAlt,iNeutral,iBlock) = &
!                      NeutralSourcesTotal(iLon,iLat,iAlt,iNeutral,iBlock) + &
!                      NeutralSources(iNeutral) * DtSub
!
!                 NeutralLossesTotal(iLon,iLat,iAlt,iNeutral,iBlock) = &
!                      NeutralLossesTotal(iLon,iLat,iAlt,iNeutral,iBlock) + &
!                      NeutralLosses(iNeutral) * DtSub


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

              EmissionTotals = EmissionTotals + Emission*DtSub

              DtTotal = DtTotal + DtSub

              if (DtSub < DtMin) DtMin = DtSub
              niters = 1


              do while (DtTotal < Dt)

                 DtSub = Dt - DtTotal

                 call calc_chemical_sources(iLon,iLat,iAlt,iBlock,IonSources,IonLosses,NeutralSources, &
                      NeutralLosses,ChemicalHeatingSub,Emission)

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

                 do iIon = 1, nIons-1
                    if (IonSources(iion) > 100.0*ions(iion) .and. IonLosses(iion) > 100.0 * ions(iion)) then
                       F = ions(iion)/IonSources(iion) * 10.0
                       IonLosses(iion) = f * IonLosses(iion)
                       IonSources(iion) = f * IonSources(iion)
                    endif
                 enddo

                 call calc_dtsub(IonSources,IonLosses,NeutralSources,NeutralLosses,dtSub)

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

!                    NeutralSourcesTotal(iLon,iLat,iAlt,iNeutral,iBlock) = &
!                         NeutralSourcesTotal(iLon,iLat,iAlt,iNeutral,iBlock) + &
!                         NeutralSources(iNeutral) * DtSub

!                    NeutralLossesTotal(iLon,iLat,iAlt,iNeutral,iBlock) = &
!                         NeutralLossesTotal(iLon,iLat,iAlt,iNeutral,iBlock) + &
!                         NeutralLosses(iNeutral) * DtSub


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

                 EmissionTotals = EmissionTotals + Emission*DtSub

                 DtTotal = DtTotal + DtSub

                 if (DtSub < DtMin) DtMin = DtSub

                 if (DtSub < 1.0e-9 .and. abs(DtTotal-Dt) > DtSub) then
                    write(*,*) "Chemistry is too fast!!", DtSub

                    ! Check Ions
                    write(*,*) "iAlt: ",iAlt," nIters: ",niters
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
              !              write(*,*)"actual: ", niters

              !               if (ialt .eq. 42 .and. istep .eq. 10000) then
              !                  write(*,*) "expl:"
              !                  do ineutral = 1, nspeciestotal
              !                     write(*,*) ialt,"end chem: ",iNeutral,neutrals(ineutral),NDensityS(1,1,ialt,ineutral,1),&
              !                          Neutrals(ineutral)-NDensityS(1,1,ialt,ineutral,1)
              !                  enddo
              !                  do iion = 1, nions - 1
              !                     write(*,*) ialt,"end chem: ",iIon,Ions(iion),IDensityS(1,1,ialt,iion,1),&
              !                          ions(iion)-IDensityS(1,1,ialt,iion,1)
              ! 
              !                  enddo
              !                  stop
              !               endif

              IDensityS(iLon,iLat,iAlt,1:nIons,iBlock) = Ions
              NDensityS(iLon,iLat,iAlt,:,iBlock) = Neutrals

              Emissions(iLon, iLat, iAlt, :, iBlock) =  &
                   Emissions(iLon, iLat, iAlt, :, iBlock) + EmissionTotals

           endif
        enddo
     enddo
  enddo
  !stop
  if (iDebugLevel > 2) &
       write(*,*) "===> calc_chemistry: Average Dt for this timestep : ", &
       (Dt*nLats*nLons*nAlts)/nIters



  if (iDebugLevel > 3) then
     do iIon = 1, nIons
        write(*,*) "====> calc_chemistry: Max Ion Density: ", iIon, &
             maxval(IDensityS(1:nLons,1:nLats,(nAlts*4)/5,iIon,iBlock))
     enddo
  endif
  call end_timing("calc_chemistry")


end subroutine calc_chemistry

subroutine calc_dtsub(IonSources,IonLosses,NeutralSources,NeutralLosses,dtSub)
  
  use ModChemistry, only: Ions,Neutrals
  use ModPlanet, only : nIons,nSpeciesTotal

  real, intent(inout) :: dtSub,IonLosses(nIons)
  real, intent(in) :: IonSources(nIons)
  real, intent(in) :: NeutralSources(nSpeciesTotal),NeutralLosses(nSpeciesTotal)

  real :: tli(nIons), tsi(nIons), tln(nSpeciesTotal), tsn(nSpeciesTotal)
  real :: ILoss(nions)
  integer :: iion

  ! Ions

  tli = DtSub * IonLosses
  tsi = DtSub * IonSources + 0.25*Ions
  
  do while (minval(tsi-tli) < 0.0 .and. DtSub > 1.0e-16)
     DtSub = DtSub/2.0
     tli = DtSub * IonLosses
     tsi = DtSub * IonSources + 0.25*Ions
  enddo
  
!  tli = DtSub * IonLosses
!  tsi = DtSub * IonSources + Ions
!
!  do iIon = 1, nIons-1
!     do while (tsi(iIon)-tli(iIon) < 0.0 .and. DtSub > 1.0e-2)
!        if (tsi(iIon)-tli(iIon) < 0.0 .and. Ions(iIon) < 1.0e7) then
!           IonLosses(iIon) = &
!                (IonSources(iIon) + Ions(iIon)/DtSub)*0.9
!        else
!           DtSub = DtSub/2.0
!        endif
!        tli(iIon) = DtSub * IonLosses(iIon)
!        tsi(iIon) = DtSub * IonSources(iIon) + Ions(iIon)
!     enddo
!  enddo

  !---- Neutrals
  tln = DtSub * NeutralLosses
  tsn = DtSub * NeutralSources + 0.1*Neutrals
  do while (minval(tsn-tln) < 0.0)
     DtSub = DtSub/2.0
     tln = DtSub * NeutralLosses
     tsn = DtSub * NeutralSources + 0.1*Neutrals
  enddo

end subroutine calc_dtsub

