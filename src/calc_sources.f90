!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine calc_GITM_sources(iBlock)

  use ModInputs
  use ModSources
  use ModGITM

  implicit none

  integer, intent(in) :: iBlock

  integer :: iAlt, iError, iDir, iLat, iLon, iSpecies
  integer :: iiAlt, iiLat, iiLon
  real :: tmp(nLons, nLats, nAlts)
  real :: tmp2(nLons, nLats, 0:nAlts+1)
  real :: tmp3(nLons, nLats, 0:nAlts+1)
  real :: RhoI(nLons, nLats, nAlts)
  real :: Rho110(nLons, nLats, 0:nAlts+1)
  real :: ScaleHeight(-1:nLons+2, -1:nLats+2, -1:nAlts+2)
  real :: Prandtl(nLons,nLats,0:nalts+1)

  ! JMB Update 07/13/2017
  ! Updated Conduction Routines Require intermediate
  ! Steps outlined in Harwood et al. [2016]
  ! Oscillation-free method for semi-linear diffusion equations
  ! digitalcommons.georgefox.edu  

  ! Temperature variables
  real :: TemperatureStage1(1:nLons,1:nLats,-1:nAlts+2)
  real :: TemperatureH(1:nLons,1:nLats,-1:nAlts+2)
  real :: TemperatureF(1:nLons,1:nLats,-1:nAlts+2)

  ! Vertical Viscosity Variables
  real :: VerticalVelocityStage0(1:nLons,1:nLats,-1:nAlts+2,1:nSpecies)
  real :: VerticalVelocityStage1(1:nLons,1:nLats,-1:nAlts+2,1:nSpecies)
  real :: VerticalVelocityH(1:nLons,1:nLats,-1:nAlts+2,1:nSpecies)
  real :: VerticalVelocityF(1:nLons,1:nLats,-1:nAlts+2,1:nSpecies)

  ! Bulk Viscosity Variables
  real :: VelocityStage0(1:nLons,1:nLats,-1:nAlts+2)
  real :: VelocityStage1(1:nLons,1:nLats,-1:nAlts+2,1:3)
  real :: VelocityH(1:nLons,1:nLats,-1:nAlts+2,1:3)
  real :: VelocityF(1:nLons,1:nLats,-1:nAlts+2,1:3)

  ! Sub-timestep used in the multi-step conduction
  real :: DtCSLocal
  ! Species-Specific Boundary Condition Settings
  logical ::VertVelNeuBCS(1:nSpecies)   ! Set to true if you want Neumann BCs
  ! NeuBCS used for Conduction and Bulk Winds
  logical ::NeuBCS                      ! Set to true if you want Neumann BCs

  ! Note:  NeuBCs is a setting for Neumann Boundary Conditions.
  ! If you are not specifying a Jeans Escape at the top, then
  ! NeuBCS = .true.
  ! 
  ! If you are setting a Jeans escape (i.e. for H, H2, or He)
  ! then NeuBCS MUST BE SET TO FALSE FOR THAT SPECIES
  ! Otherwise, you will not get the escape that you want!

  call report("calc_GITM_sources",1)

  ! calc_rate is used to determine reaction rates, heating coefficients,
  ! ion-neutral collision frequency, lambdas for ion drag, etc.

  if (iDebugLevel > 4) write(*,*) "=====> going into calc_rates", iproc
  ! JMB: 07/13/2017
  ! Initialize our 2nd Order Multi-Stage Variables
  TemperatureStage1(1:nLons,1:nLats,-1:nAlts+2) = &
       Temperature(1:nLons,1:nLats,-1:nAlts+2,iBlock)
  TemperatureH(1:nLons,1:nLats,-1:nAlts+2) = &
       Temperature(1:nLons,1:nLats,-1:nAlts+2,iBlock)
  TemperatureF(1:nLons,1:nLats,-1:nAlts+2) = &
       Temperature(1:nLons,1:nLats,-1:nAlts+2,iBlock)

  VerticalVelocityStage1(1:nLons,1:nLats,-1:nAlts+2,1:nSpecies) = &
       VerticalVelocity(1:nLons,1:nLats,-1:nAlts+2,1:nSpecies,iBlock)
  VerticalVelocityH(1:nLons,1:nLats,-1:nAlts+2,1:nSpecies) = &
       VerticalVelocity(1:nLons,1:nLats,-1:nAlts+2,1:nSpecies,iBlock)
  VerticalVelocityF(1:nLons,1:nLats,-1:nAlts+2,1:nSpecies) = &
       VerticalVelocity(1:nLons,1:nLats,-1:nAlts+2,1:nSpecies,iBlock)

  VelocityStage1(1:nLons,1:nLats,-1:nAlts+2,1:3) = &
       Velocity(1:nLons,1:nLats,-1:nAlts+2,1:3,iBlock)
  VelocityH(1:nLons,1:nLats,-1:nAlts+2,1:3) = &
       Velocity(1:nLons,1:nLats,-1:nAlts+2,1:3,iBlock)
  VelocityF(1:nLons,1:nLats,-1:nAlts+2,1:3) = &
       Velocity(1:nLons,1:nLats,-1:nAlts+2,1:3,iBlock)

  ! Note, this assumes that all species need Neumann Upper boundary
  ! conditions.  If you specify escape rates in set_vertical_bcs.f90
  ! such as Jeans Escape, then set VertVelNeuBCs(iH_ or iH2_) = .false.
  ! This will allow you to specify escape rates at the top.
  VertVelNeuBCS(1:nSpecies) = .true. 

  ChemicalHeatingRate = 0.0
  ChemicalHeatingRateIon = 0.0
  ChemicalHeatingRateEle = 0.0

  ChemicalHeatingSpecies = 0.0

  call calc_eddy_diffusion_coefficient(iBlock)  
  call calc_rates(iBlock)
  call calc_collisions(iBlock)

  RhoI = IDensityS(1:nLons,1:nLats,1:nAlts,ie_,iBlock) * &
       MeanIonMass(1:nLons,1:nLats,1:nAlts)

  !\
  ! ---------------------------------------------------------------
  ! These terms are for Neutral Temperature
  ! ---------------------------------------------------------------
  !/

  !\
  ! Solar Heating -------------------------------------------------
  !/

  if (iDebugLevel > 4) write(*,*) "=====> solar heating", iproc
  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)

  if (UseSolarHeating .or. UseIonChemistry) then

     ! So far, calc_physics only has stuff that is needed for solar
     ! euv, such as solar zenith angles, and local time.

     call calc_physics(iBlock)
     call euv_ionization_heat(iBlock)

  endif

  if (.not. UseSolarHeating) EuvHeating = 0.0

  ! This includes Radiative Cooling....
  RadCooling = 0.0
  call calc_planet_sources(iBlock)

  !\
  ! Conduction ----------------------------------------------------
  !/

  if (iDebugLevel > 4) write(*,*) "=====> conduction", iproc
  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)

  Rho110 = Rho(1:nLons, 1:nLats,0:nAlts+1, iBlock)
  
  if(UseConduction)then

     ! Note: Neumann = .true. is needed if you want 
     ! the gradient to be zero at the top.

     NeuBCS = .true.  ! Use Neumann Boundary Conditions

     tmp2 = Rho110 * cp(1:nLons, 1:nLats,0:nAlts+1, iBlock)

     if (UseTurbulentCond) then
        Prandtl = &
             KappaEddyDiffusion(1:nLons, 1:nLats,0:nAlts+1, iBlock) * &
                            Rho110 * &
                             Cp(1:nLons, 1:nLats,0:nAlts+1, iBlock)
     else 
        Prandtl = 0.0
     endif

     DtCSLocal = Dt/2.0
     call calc_conduction(iBlock, DtCSLocal, NeuBCS, &
          Temperature(1:nLons, 1:nLats,-1:nAlts+2, iBlock) * &
             TempUnit(1:nLons, 1:nLats,-1:nAlts+2), &
            KappaTemp(1:nLons, 1:nLats, 0:nAlts+1, iBlock) + &
              Prandtl(1:nLons, 1:nLats, 0:nAlts+1), &
                 tmp2(1:nLons, 1:nLats, 0:nAlts+1), &
       MoleConduction(1:nLons, 1:nLats, 0:nAlts+1))

     TemperatureStage1(1:nLons,1:nLats,0:nAlts+1) = &
           Temperature(1:nLons,1:nLats,0:nAlts+1,iBlock) + &
        MoleConduction(1:nLons,1:nLats,0:nAlts+1)/&
              TempUnit(1:nLons,1:nLats,0:nAlts+1)

     DtCSLocal = Dt/2.0
     call calc_conduction(iBlock, DtCSLocal, NeuBCS, &
          TemperatureStage1(1:nLons,1:nLats,-1:nAlts+2) * &
                   TempUnit(1:nLons,1:nLats,-1:nAlts+2), &
                  KappaTemp(1:nLons,1:nLats,0:nAlts+1,iBlock) + &
                    Prandtl(1:nLons,1:nLats,0:nAlts+1), &
                       tmp2(1:nLons,1:nLats,0:nAlts+1), &
             MoleConduction(1:nLons,1:nLats,0:nAlts+1))

     TemperatureH(1:nLons,1:nLats, 0:nAlts+1) = &
          TemperatureStage1(1:nLons,1:nLats,0:nAlts+1) + &
             MoleConduction(1:nLons,1:nLats,0:nAlts+1)/&
                   TempUnit(1:nLons,1:nLats,0:nAlts+1)

     ! Full Time Step Update
     DtCSLocal = Dt
     call calc_conduction(iBlock, DtCSLocal, NeuBCS, &
          Temperature(1:nLons, 1:nLats,-1:nAlts+2, iBlock) * &
             TempUnit(1:nLons, 1:nLats,-1:nAlts+2), &
            KappaTemp(1:nLons, 1:nLats, 0:nAlts+1, iBlock) + &
              Prandtl(1:nLons, 1:nLats, 0:nAlts+1 ), &
                 tmp2(1:nLons, 1:nLats, 0:nAlts+1 ), &
       MoleConduction(1:nLons, 1:nLats, 0:nAlts+1))

     TemperatureF(1:nLons,1:nLats, 0:nAlts+1) = &
          Temperature(1:nLons,1:nLats, 0:nAlts+1,iBlock) + &
       MoleConduction(1:nLons,1:nLats, 0:nAlts+1) / &
             TempUnit(1:nLons,1:nLats, 0:nAlts+1)

     ! Note that Conduction is the net temperature update
     Conduction(1:nLons,1:nLats,0:nAlts+1) = &
          (2.0*TemperatureH(1:nLons,1:nLats,0:nAlts+1) - &
               TemperatureF(1:nLons,1:nLats,0:nAlts+1)) - &
                Temperature(1:nLons,1:nLats,0:nAlts+1,iBlock) 

  else
     Conduction = 0.0
  end if

  !\
  ! ---------------------------------------------------------------
  ! These terms are for Neutral Winds
  ! ---------------------------------------------------------------
  !/

  if (iDebugLevel > 4) write(*,*) "=====> Ion Drag", iproc
  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  
  if (UseIonDrag) then

     tmp = Collisions(1:nLons,1:nLats,1:nAlts,iVIN_)*&
          RhoI/Rho(1:nLons,1:nLats,1:nAlts,iBlock)

     do iDir = 1, 3
        IonDrag(:,:,:,iDir) = tmp * &
             (IVelocity(1:nLons,1:nLats,1:nAlts,iDir,iBlock) - &
             Velocity(1:nLons,1:nLats,1:nAlts,iDir,iBlock))
     enddo

     ! F_iondrag = rho_i/rho * Vis * (Ui-Un)
     ! where Vis = Vin *(Ns/N)  

     do iSpecies = 1, nSpecies

        tmp = Collisions(1:nLons,1:nLats,1:nAlts,iVIN_)*&
             RhoI / &
             (Mass(iSpecies) * &
             NDensityS(1:nLons,1:nLats,1:nAlts,iSpecies,iBlock)) * &
             (NDensityS(1:nLons,1:nLats,1:nAlts,iSpecies,iBlock) / &
             NDensity(1:nLons,1:nLats,1:nAlts,iBlock))

        VerticalIonDrag(:,:,:,iSpecies) = tmp * &
             (IVelocity(1:nLons,1:nLats,1:nAlts,iUp_,iBlock) - &
             VerticalVelocity(1:nLons,1:nLats,1:nAlts,iSpecies,iBlock))

     enddo

  else

     IonDrag = 0.0
     VerticalIonDrag = 0.0

  endif

  !\
  ! Viscosity ----------------------------------------------------
  !/

  if (iDebugLevel > 4) write(*,*) "=====> Viscosity", iproc
  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  
  if (UseViscosity) then

     if (iDebugLevel > 4) write(*,*) "=====> calc_viscosity", iproc
     if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
     call calc_viscosity(iBlock)

     ! 07/13/2017:  JMB, Below is the 2nd Order
     ! Implicit method for Horizontal Winds.
     ! We assume that these winds (North, East)
     ! have a zero gradient at the top, so NeuBCS = .true.

     NeuBCS = .true.
     do iDir = iEast_, iNorth_

        if (iDebugLevel > 4) write(*,*) "=====> calc_conduct", iproc, iDir, 1
        if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)

        DtCSLocal = Dt/2.0
        VelocityStage0 = Velocity(1:nLons, 1:nLats,-1:nAlts+2, iDir,iBlock)
        call calc_conduction(iBlock, DtCSLocal, NeuBCS, &
            VelocityStage0, & 
            ViscCoef(1:nLons, 1:nLats, 0:nAlts+1      ) + &
                 Rho110 * &
  KappaEddyDiffusion(1:nLons, 1:nLats, 0:nAlts+1, iBlock) , &
                 Rho110, &
           Viscosity(1:nLons, 1:nLats, 0:nAlts+1, iDir))

     VelocityStage1(1:nLons,1:nLats, 0:nAlts+1, iDir ) = &
           Velocity(1:nLons,1:nLats, 0:nAlts+1, iDir ,iBlock) + &
          Viscosity(1:nLons,1:nLats, 0:nAlts+1, iDir  )

        if (iDebugLevel > 4) write(*,*) "=====> calc_conduct", iproc, iDir, 2
        if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)

        DtCSLocal = Dt/2.0
        call calc_conduction(iBlock, DtCSLocal, NeuBCS, &
              VelocityStage1(1:nLons, 1:nLats,-1:nAlts+2, iDir), & 
                    ViscCoef(1:nLons, 1:nLats, 0:nAlts+1      ) + &
                         Rho110 * &
          KappaEddyDiffusion(1:nLons, 1:nLats, 0:nAlts+1, iBlock) , &
                         Rho110, &
                   Viscosity(1:nLons, 1:nLats, 0:nAlts+1,iDir))

          VelocityH(1:nLons,1:nLats, 0:nAlts+1, iDir ) = &
     VelocityStage1(1:nLons,1:nLats, 0:nAlts+1, iDir ) + &
          Viscosity(1:nLons,1:nLats, 0:nAlts+1, iDir  )

     if (iDebugLevel > 4) write(*,*) "=====> calc_conduct", iproc, iDir, 3
     if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)

     DtCSLocal = Dt
     VelocityStage0 = Velocity(1:nLons, 1:nLats,-1:nAlts+2, iDir,iBlock)
     call calc_conduction(iBlock, DtCSLocal, NeuBCS, &
            VelocityStage0, & 
            ViscCoef(1:nLons, 1:nLats, 0:nAlts+1      ) + &
                 Rho110 * &
  KappaEddyDiffusion(1:nLons, 1:nLats, 0:nAlts+1, iBlock) , &
                 Rho110, &
           Viscosity(1:nLons, 1:nLats, 0:nAlts+1,iDir))

     VelocityF(1:nLons,1:nLats, 0:nAlts+1, iDir) = &
      Velocity(1:nLons,1:nLats, 0:nAlts+1, iDir,iBlock) + &
     Viscosity(1:nLons,1:nLats, 0:nAlts+1, iDir  )

      ! This is the final 2nd ORder Implicit Viscosity for 
      ! the horizontal Bulk Winds
        Viscosity(1:nLons,1:nLats, 0:nAlts+1, iDir ) = &
  ( 2.0*VelocityH(1:nLons,1:nLats, 0:nAlts+1, iDir         ) - &
        VelocityF(1:nLons,1:nLats, 0:nAlts+1, iDir         )) - &
         Velocity(1:nLons,1:nLats, 0:nAlts+1,iDir,iBlock) 

     enddo !iDir = iEast_, iNorth_

     if (iDebugLevel > 4) write(*,*) "=====> Vertical Viscosity", iproc
     if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)

     ! JMB:  07/13/2017
     ! Add 2nd Order Viscosity for Each Species
     Viscosity(:,:,:,iUp_) = 0.0
     do iSpecies = 1, nSpecies

          DtCSLocal = Dt/2.0
          VelocityStage0 = VerticalVelocity(1:nLons, 1:nLats,-1:nAlts+2, iSpecies,iBlock)
          call calc_conduction(iBlock, DtCSLocal, VertVelNeuBCS(iSpecies), &
               VelocityStage0, & 
            (4.0/3.0)*ViscCoefS(1:nLons, 1:nLats, 0:nAlts+1, iSpecies) + &
       Mass(iSpecies)*NDensityS(1:nLons, 1:nLats, 0:nAlts+1, iSpecies, iBlock) * &
             KappaEddyDiffusion(1:nLons, 1:nLats, 0:nAlts+1, iBlock) , &
       Mass(iSpecies)*NDensityS(1:nLons, 1:nLats, 0:nAlts+1, iSpecies, iBlock), &
                      Viscosity(1:nLons, 1:nLats, 0:nAlts+1, iUp_))

          VerticalVelocityStage1(1:nLons,1:nLats, 0:nAlts+1, iSpecies ) = &
                VerticalVelocity(1:nLons,1:nLats, 0:nAlts+1, iSpecies ,iBlock) + &
                       Viscosity(1:nLons,1:nLats, 0:nAlts+1, iUp_  )

          DtCSLocal = Dt/2.0
          call calc_conduction(iBlock, DtCSLocal, VertVelNeuBCS(iSpecies), &
         VerticalVelocityStage1(1:nLons, 1:nLats,-1:nAlts+2, iSpecies), & 
            (4.0/3.0)*ViscCoefS(1:nLons, 1:nLats, 0:nAlts+1, iSpecies) + &
       Mass(iSpecies)*NDensityS(1:nLons, 1:nLats, 0:nAlts+1, iSpecies, iBlock) * &
             KappaEddyDiffusion(1:nLons, 1:nLats, 0:nAlts+1, iBlock) , &
       Mass(iSpecies)*NDensityS(1:nLons, 1:nLats, 0:nAlts+1, iSpecies, iBlock), &
                      Viscosity(1:nLons, 1:nLats, 0:nAlts+1, iUp_))

               VerticalVelocityH(1:nLons,1:nLats, 0:nAlts+1, iSpecies ) = &
          VerticalVelocityStage1(1:nLons,1:nLats, 0:nAlts+1, iSpecies ) + &
                       Viscosity(1:nLons,1:nLats, 0:nAlts+1, iUp_  )

          ! Full Time Step Update
          DtCSLocal = Dt
          VelocityStage0 = VerticalVelocity(1:nLons, 1:nLats,-1:nAlts+2, iSpecies,iBlock)
          call calc_conduction(iBlock, DtCSLocal, VertVelNeuBCS(iSpecies), &
               VelocityStage0, &
            (4.0/3.0)*ViscCoefS(1:nLons, 1:nLats, 0:nAlts+1, iSpecies) + &
       Mass(iSpecies)*NDensityS(1:nLons, 1:nLats, 0:nAlts+1, iSpecies, iBlock) * &
             KappaEddyDiffusion(1:nLons, 1:nLats, 0:nAlts+1, iBlock) , &
       Mass(iSpecies)*NDensityS(1:nLons, 1:nLats, 0:nAlts+1, iSpecies, iBlock), &
                      Viscosity(1:nLons, 1:nLats, 0:nAlts+1, iUp_))

          VerticalVelocityF(1:nLons,1:nLats, 0:nAlts+1, iSpecies) = &
           VerticalVelocity(1:nLons,1:nLats, 0:nAlts+1, iSpecies,iBlock) + &
                  Viscosity(1:nLons,1:nLats, 0:nAlts+1, iUp_  )

          ! Simply update the Vertical Winds here (rather than in add_sources)
               VerticalVelocity(1:nLons,1:nLats, 0:nAlts+1, iSpecies, iBlock  ) = &
          2.0*VerticalVelocityH(1:nLons,1:nLats, 0:nAlts+1, iSpecies         ) - &
              VerticalVelocityF(1:nLons,1:nLats, 0:nAlts+1, iSpecies         ) 

      ! This is the final 2nd ORder Implicit Viscosity for 
      ! the horizontal Bulk Winds
        VerticalViscosityS(1:nLons,1:nLats, 0:nAlts+1, iSpecies ) = &
  ( 2.0*VelocityH(1:nLons,1:nLats, 0:nAlts+1, iDir         ) - &
        VelocityF(1:nLons,1:nLats, 0:nAlts+1, iDir         )) - &
         VerticalVelocity(1:nLons,1:nLats, 0:nAlts+1,iSpecies,iBlock) 

     enddo ! iSpecies

  else
     Viscosity = 0.0
  end if
!
  !\
  ! ---------------------------------------------------------------
  ! These terms are for Ion Densities
  ! ---------------------------------------------------------------
  !/

  if (iDebugLevel > 4) write(*,*) "=====> get_potential", iproc
  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  call get_potential(iBlock)

  if (iDebugLevel > 4) write(*,*) "=====> Efield", iproc
  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  call calc_efield(iBlock)

  if (iDebugLevel > 4) write(*,*) "=====> Aurora", iproc
  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  call aurora(iBlock)
  if (UseAuroralHeating) then
     AuroralHeating = AuroralHeatingRate(:,:,:,iBlock) / &
          TempUnit(1:nLons,1:nLats,1:nAlts) / cp(:,:,1:nAlts,iBlock) / &
          rho(1:nLons,1:nLats,1:nAlts, iBlock)
  else
     AuroralHeating = 0.0
  endif

  if (iDebugLevel > 4) write(*,*) "=====> Ion Velocity", iproc
  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  call calc_ion_v(iBlock)

  if (iDebugLevel > 4) write(*,*) "=====> Chemistry", iproc
  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  ! The Emissions array was never set. Should this be here or earlier ????
  Emissions(:,:,:,:,iBlock) = 0.0
  call calc_chemistry(iBlock)

  ChemicalHeatingRate(:,:,:) = &
       ChemicalHeatingRate(:,:,:) * Element_Charge / &
       TempUnit(1:nLons,1:nLats,1:nAlts) / cp(1:nLons,1:nLats,1:nAlts,iBlock)/&
       rho(1:nLons,1:nLats,1:nAlts,iBlock)
	   
  ChemicalHeatingRateIon(:,:,:) = &
       ChemicalHeatingRateIon(:,:,:) * Element_Charge

  ChemicalHeatingSpecies = ChemicalHeatingSpecies * Element_Charge

end subroutine calc_GITM_sources
