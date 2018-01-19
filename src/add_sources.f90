!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine add_sources

  use ModGITM
  use ModTime
  use ModSources
  use ModInputs
  use ModUserGITM

  implicit none

  integer :: iBlock, iLon, iLat, iAlt, iSpecies
  integer :: iDir
  logical :: IsFirstTime=.true.

  real(kind=8), dimension(0:nLons+1,0:nLats+1,0:nAlts+1) :: &
       eHeatingp, iHeatingp, eHeatingm, iHeatingm, iHeating, lame, lami

  call report("add_sources",2)

  if (floor((tSimulation-dt)/DtPotential) /= &
       floor((tsimulation)/DtPotential) .or. IsFirstTime) then
     if (UseDynamo .and. .not. Is1D) then
        call UA_calc_electrodynamics(iLon, iLat)
     else
        call UA_calc_electrodynamics_1d
     endif
     IsFirstTime = .false.
  endif

  do iBlock = 1, nBlocks

     ! All the physics is left out or added in in calc_GITM_sources.  If
     ! you want to turn something off, look for the UseWhatever variable
     ! in calc_GITM_sources.  Then fill the source with 0.0, so this routine
     ! does not change.

     call calc_GITM_sources(iBlock)
     call calc_electron_ion_sources(iBlock, &
          eHeatingp,iHeatingp,eHeatingm,iHeatingm,iHeating,lame,lami)

     !! To turn off EuvHeating, turn UseSolarHeating=.false. in UAM.in
     !! To turn off JouleHeating, turn UseJouleHeating=.false. in UAM.in
     !! To turn off AuroralHeating, turn Use=AuroralHeating.false. in UAM.in
     !! To turn off Conduction, turn UseConduction=.false. in UAM.in

     ! JMB:  07/13/2017.
     ! 2nd order conduction update:  Separately add Conduction to this
     ! because Conduction now spans(0:nAlts+1)
     Temperature(1:nLons, 1:nLats, 1:nAlts, iBlock) = &
          Temperature(1:nLons, 1:nLats, 1:nAlts, iBlock) + Dt * ( &
          LowAtmosRadRate(1:nLons, 1:nLats, 1:nAlts, iBlock) &
          /TempUnit(1:nLons,1:nLats,1:nAlts)&
          - RadCooling(1:nLons, 1:nLats, 1:nAlts, iBlock) &
          + EuvHeating(1:nLons, 1:nLats, 1:nAlts, iBlock) &
          + PhotoElectronHeating(1:nLons, 1:nLats, 1:nAlts, iBlock) &
          + AuroralHeating &
          + JouleHeating &
          + ElectronHeating &
          ) &
          + ChemicalHeatingRate &
          + UserHeatingRate(1:nLons, 1:nLats, 1:nAlts, iBlock)

     Temperature(1:nLons, 1:nLats, 0:nAlts+1, iBlock) = &
          Temperature(1:nLons, 1:nLats, 0:nAlts+1, iBlock) + &
          + Conduction(1:nLons,1:nLats,0:nAlts+1)

     !-------------------------------------------
     ! This is an example of a user output:
 
     UserData3D(:,:,:,1,iBlock) = 0.0
     UserData3D(1:nLons, 1:nLats, 1:nAlts, 1, iBlock) = JouleHeating

     do while (minval(temperature(1:nLons, 1:nLats, 1:nAlts, iBlock)) < 0.0)
        write(*,*) "Negative Temperature Found!!!  Correcting!!!"
        do iLon = 1, nLons
           do iLat = 1, nLats
              iAlt = 1
              if (temperature(iLon, iLat, iAlt, iBlock) < 0.0) &
                   temperature(iLon, iLat, iAlt, iBlock) = &
                   temperature(iLon, iLat, iAlt-1, iBlock)
              do iAlt = 2, nAlts
                 if (temperature(iLon, iLat, iAlt, iBlock) < 0.0) then
                    temperature(iLon, iLat, iAlt, iBlock) = &
                         (temperature(iLon, iLat, iAlt-1, iBlock) +  &
                         temperature(iLon, iLat, iAlt+1, iBlock))/2.0

                    write(*,*) "Sources : ", &
                         temperature(iLon, iLat, iAlt, iBlock), &
                         EuvHeating(iLon, iLat, iAlt, iBlock) * dt, &
                         RadCooling(iLon, iLat, iAlt, iBlock) * dt, &
                         AuroralHeating(iLon, iLat, iAlt) * dt, &
                         JouleHeating(iLon, iLat, iAlt) * dt, &
                         Conduction(iLon, iLat, iAlt), &
                         ChemicalHeatingRate(iLon, iLat, iAlt)
                    call stop_gitm('Negative Temperature Found')

                 endif
              enddo
           enddo
        enddo
     enddo

     if (iDebugLevel > 2 .and. Is1D) then
!        do iAlt = 1,nAlts
iAlt = 10
           write(*,*) "===> Temp Sources : ", iAlt, dt,&
                EuvHeating(1, 1, iAlt, iBlock)*dt, &
!                NOCooling(1,1,iAlt)*dt, &
!                OCooling(1,1,iAlt)*dt, &
                AuroralHeating(1,1,iAlt)*dt, &
                JouleHeating(1,1,iAlt)*dt, &
                ChemicalHeatingRate(1,1,iAlt), &
                Conduction(1,1,iAlt), temperature(1,1,iAlt,iBlock)
!        enddo
     endif

     iAlt = nAlts-2
     if (iDebugLevel > 2) &
          write(*,*) "===> Sum Temp Sources : ", &
          sum(EuvHeating(1:nLons, 1:nLats, iAlt, iBlock))*dt, &
          sum(NOCooling(:,:,iAlt))*dt, &
          sum(OCooling(:,:,iAlt))*dt, &
          sum(AuroralHeating(:,:,iAlt))*dt, &
          sum(JouleHeating(:,:,iAlt))*dt, &
          sum(Conduction(:,:,iAlt))

     !! To turn off IonDrag, turn UseIonDrag=.false. in UAM.in
     do iDir = 1, 3
       Velocity(1:nLons, 1:nLats, 1:nAlts, iDir, iBlock) = &
            Velocity(1:nLons, 1:nLats, 1:nAlts, iDir, iBlock) + &
            Dt*IonDrag(:,:,:,iDir) + GWAccel(:,:,:,iDir)

       Velocity(1:nLons, 1:nLats, 0:nAlts+1, iDir, iBlock) = &
            Velocity(1:nLons, 1:nLats, 0:nAlts+1, iDir, iBlock) + &
            Viscosity(1:nLons,1:nLats, 0:nAlts+1,iDir)
     enddo 

     !! To turn off IonDrag, turn UseIonDrag=.false. in UAM.in
     !! To turn off NeutralFriction, turn UseNeutralFriction=.false. in UAM.in

     do iSpecies = 1, nSpecies
        VerticalVelocity(1:nLons, 1:nLats, 1:nAlts, iSpecies, iBlock) =&
             VerticalVelocity(1:nLons, 1:nLats, 1:nAlts, iSpecies, iBlock) + &
             Dt*(VerticalIonDrag(:,:,:,iSpecies)) + &
             NeutralFriction(:,:,:,iSpecies) 

        VerticalVelocity(1:nLons, 1:nLats, 0:nAlts+1, iSpecies, iBlock) =&
        VerticalVelocity(1:nLons, 1:nLats, 0:nAlts+1, iSpecies, iBlock) +&
                VerticalViscosityS(1:nLons,1:nLats,0:nAlts+1,iSpecies)

     enddo

    
     call calc_electron_temperature(iBlock,eHeatingp,iHeatingp,eHeatingm,iHeatingm,iHeating,lame,lami)

     !! To turn off Diffusion, turn UseDiffusion=.false. in UAM.in
     do iLon = 1, nLons
        do iLat = 1, nLats
           do iAlt = 1, nAlts
              Rho(iLon, iLat, iAlt, iBlock) = &
                   sum(Mass(1:nSpecies) * &
                   NDensityS(iLon,iLat,iAlt,1:nSpecies,iBlock) )
              NDensity(iLon, iLat, iAlt, iBlock) = &
                   sum(NDensityS(iLon,iLat,iAlt,1:nSpecies,iBlock) )
           enddo
        enddo
     enddo

     Velocity(1:nLons, 1:nLats, 1:nAlts, iUp_, iBlock) = 0.0
     do iSpecies = 1, nSpecies
        Velocity(1:nLons, 1:nLats, 1:nAlts, iUp_, iBlock) = &
             Velocity(1:nLons, 1:nLats, 1:nAlts, iUp_, iBlock) + &
             VerticalVelocity(1:nLons, 1:nLats, 1:nAlts, iSpecies, iBlock)* &
             Mass(iSpecies) * &
             NDensityS(1:nLons, 1:nLats, 1:nAlts, iSpecies, iBlock) / &
             Rho(1:nLons, 1:nLats, 1:nAlts, iBlock)
     enddo

  enddo

end subroutine add_sources
