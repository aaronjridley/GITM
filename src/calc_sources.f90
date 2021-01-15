!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine calc_GITM_sources(iBlock)

  use ModInputs
  use ModSources
  use ModGITM
  use ModEUV, only: sza
  use ModTime, only: iTimeArray, UTime
  use ModPlanet, only: IsVenus

  implicit none

  integer, intent(in) :: iBlock

  integer :: iAlt, iError, iDir, iLat, iLon, iSpecies, iIon
  integer :: iiAlt, iiLat, iiLon, iSza, jAlt
  real :: m !slope for IR heating
  real :: tmp(nLons, nLats, nAlts)
  real :: tmp2(nLons, nLats, 0:nAlts+1)
  real :: tmp3(nLons, nLats, 0:nAlts+1)
  real :: RhoI(nLons, nLats, nAlts)
  real :: Rho110(nLons, nLats, 0:nAlts+1)
  real :: Weight(nLons, nLats, nAlts)
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

  real :: TmpTemp(nLons, nLats, -1:nAlts+2)
  real :: TmpDiff(nLons, nLats, 0:nAlts+1)
  real :: TmpMulFac(nLons, nLats, 0:nAlts+1)

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

  real :: TmpVel(nLons, nLats, -1:nAlts+2)
  real :: CondResult(nLons, nLats, 0:nAlts+1)

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

  !BP: Proportionality for IR interpolation in sza and in altitude
  real :: r_sza
  real :: x03, x12
  real :: LST
  real :: RP_hours = Rotation_Period/3600.0
  integer :: iLatTable, iAltTable

  !BP: sza data for IR table                                                                 
  real, dimension(11) :: sza_table = (/0.0,25.0,36.0,45.0,53.0,60.0, &
                                      66.0,72.0,78.0,84.0,89.0 /)*pi/180.0
  real, dimension(40) :: altitude_table = (/80.0,82.05,84.10,86.15,88.21, &
                                            90.26,92.31,94.36,96.41,98.46, &
                                            100.51,102.56,104.62,106.67,108.72, &
                                            110.77,112.82,114.87,116.92,118.97, &
                                            121.03,123.08,125.13,127.18,129.23, &
                                            131.28,133.33,135.38,137.44,139.49, &
                                            141.54,143.59,145.64,147.69,149.74, &
                                            151.79,153.85,155.90,157.95,160.00/)*1000.0


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
     TmpTemp = Temperature(1:nLons, 1:nLats,-1:nAlts+2, iBlock) * &
          TempUnit(1:nLons, 1:nLats,-1:nAlts+2)
     TmpDiff = KappaTemp(1:nLons, 1:nLats, 0:nAlts+1, iBlock) + &
          Prandtl(1:nLons, 1:nLats, 0:nAlts+1)

     call calc_conduction(iBlock, DtCSLocal, NeuBCS, &
          TmpTemp, TmpDiff, tmp2, MoleConduction)

     TemperatureStage1(1:nLons,1:nLats,0:nAlts+1) = &
           Temperature(1:nLons,1:nLats,0:nAlts+1,iBlock) + &
        MoleConduction(1:nLons,1:nLats,0:nAlts+1)/&
              TempUnit(1:nLons,1:nLats,0:nAlts+1)

     DtCSLocal = Dt/2.0
     TmpTemp = TemperatureStage1(1:nLons, 1:nLats,-1:nAlts+2) * &
          TempUnit(1:nLons, 1:nLats,-1:nAlts+2)
     TmpDiff = KappaTemp(1:nLons, 1:nLats, 0:nAlts+1, iBlock) + &
          Prandtl(1:nLons, 1:nLats, 0:nAlts+1)

     call calc_conduction(iBlock, DtCSLocal, NeuBCS, &
          TmpTemp, TmpDiff, tmp2, MoleConduction)

     TemperatureH(1:nLons,1:nLats, 0:nAlts+1) = &
          TemperatureStage1(1:nLons,1:nLats,0:nAlts+1) + &
             MoleConduction(1:nLons,1:nLats,0:nAlts+1)/&
                   TempUnit(1:nLons,1:nLats,0:nAlts+1)

     ! Full Time Step Update
     DtCSLocal = Dt
     TmpTemp = Temperature(1:nLons, 1:nLats,-1:nAlts+2, iBlock) * &
          TempUnit(1:nLons, 1:nLats,-1:nAlts+2)
     TmpDiff = KappaTemp(1:nLons, 1:nLats, 0:nAlts+1, iBlock) + &
          Prandtl(1:nLons, 1:nLats, 0:nAlts+1)

     call calc_conduction(iBlock, DtCSLocal, NeuBCS, &
          TmpTemp, TmpDiff, tmp2, MoleConduction)

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
  
  IonDrag = 0.0
  VerticalIonDrag = 0.0

  !
  ! This method should work, but seems to be way, way too strong, and I don't know why:
  !

!  if (UseIonDrag) then
!
!     do iSpecies = 1, nSpecies
!
!        Weight = Mass(iSpecies) * NDensityS(1:nLons,1:nLats,1:nAlts,iSpecies,iBlock) / &
!             Rho(1:nLons,1:nLats,1:nAlts,iBlock)
!
!        do iIon = 1, nIons-1
!
!           ! Generalize to use all collisions, AJR - Nov. 3, 2019
!
!           tmp = &
!                IonCollisions(1:nLons,1:nLats,1:nAlts,iIon,iSpecies) * &
!                IDensityS(1:nLons,1:nLats,1:nAlts,iIon,iBlock) * MassI(iIon) / &
!                (NDensityS(1:nLons,1:nLats,1:nAlts,iSpecies,iBlock) * Mass(iSpecies))
!
!           do iDir = 1, 3
!              IonDrag(:,:,:,iDir) = IonDrag(:,:,:,iDir) + &
!                   tmp * &
!                   (IVelocity(1:nLons,1:nLats,1:nAlts,iDir,iBlock) - &
!                   Velocity(1:nLons,1:nLats,1:nAlts,iDir,iBlock)) * &
!                   Weight
!           enddo
!
!           VerticalIonDrag(:,:,:,iSpecies) = VerticalIonDrag(:,:,:,iSpecies) + &
!                tmp * &
!                (IVelocity(1:nLons,1:nLats,1:nAlts,iUp_,iBlock) - &
!                VerticalVelocity(1:nLons,1:nLats,1:nAlts,iSpecies,iBlock))
!
!        enddo
!     enddo
!
!  endif

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
        TmpDiff = ViscCoef(1:nLons, 1:nLats, 0:nAlts+1) + &
             Rho110 * KappaEddyDiffusion(1:nLons, 1:nLats, 0:nAlts+1, iBlock)

        call calc_conduction(iBlock, DtCSLocal, NeuBCS, &
             VelocityStage0, TmpDiff, Rho110, CondResult)

        Viscosity(1:nLons, 1:nLats, 0:nAlts+1, iDir) = CondResult

        VelocityStage1(1:nLons,1:nLats, 0:nAlts+1, iDir ) = &
             Velocity(1:nLons,1:nLats, 0:nAlts+1, iDir ,iBlock) + &
             Viscosity(1:nLons,1:nLats, 0:nAlts+1, iDir  )

        if (iDebugLevel > 4) write(*,*) "=====> calc_conduct", iproc, iDir, 2
        if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)

        DtCSLocal = Dt/2.0
        TmpVel = VelocityStage1(1:nLons, 1:nLats,-1:nAlts+2, iDir)
        TmpDiff = ViscCoef(1:nLons, 1:nLats, 0:nAlts+1) + &
             Rho110 * KappaEddyDiffusion(1:nLons, 1:nLats, 0:nAlts+1, iBlock)

        call calc_conduction(iBlock, DtCSLocal, NeuBCS, &
             TmpVel, TmpDiff, Rho110, CondResult)

        Viscosity(1:nLons, 1:nLats, 0:nAlts+1,iDir) = CondResult

        VelocityH(1:nLons,1:nLats, 0:nAlts+1, iDir ) = &
             VelocityStage1(1:nLons,1:nLats, 0:nAlts+1, iDir ) + &
             Viscosity(1:nLons,1:nLats, 0:nAlts+1, iDir  )

        if (iDebugLevel > 4) write(*,*) "=====> calc_conduct", iproc, iDir, 3
        if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)

        DtCSLocal = Dt
        VelocityStage0 = Velocity(1:nLons, 1:nLats,-1:nAlts+2, iDir,iBlock)
        TmpDiff = ViscCoef(1:nLons, 1:nLats, 0:nAlts+1) + &
             Rho110 * KappaEddyDiffusion(1:nLons, 1:nLats, 0:nAlts+1, iBlock)
        call calc_conduction(iBlock, DtCSLocal, NeuBCS, &
             VelocityStage0, TmpDiff, Rho110, CondResult)

        Viscosity(1:nLons, 1:nLats, 0:nAlts+1,iDir) = CondResult

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


     ! I think that the code below is not correct for the viscosity in the vertical winds.
     ! I need to write a viscosity solver for the vertical winds in the horizontal direction.
     
     do iSpecies = 1, nSpecies

        DtCSLocal = Dt/2.0
        VelocityStage0 = VerticalVelocity(1:nLons, 1:nLats,-1:nAlts+2, iSpecies,iBlock)
        TmpDiff = (4.0/3.0)*ViscCoefS(1:nLons, 1:nLats, 0:nAlts+1, iSpecies) + &
             Mass(iSpecies)*NDensityS(1:nLons, 1:nLats, 0:nAlts+1, iSpecies, iBlock) * &
             KappaEddyDiffusion(1:nLons, 1:nLats, 0:nAlts+1, iBlock)
        TmpMulFac = Mass(iSpecies)*NDensityS(1:nLons, 1:nLats, 0:nAlts+1, iSpecies, iBlock)

        call calc_conduction(iBlock, DtCSLocal, VertVelNeuBCS(iSpecies), &
             VelocityStage0, TmpDiff, TmpMulFac, CondResult)

        Viscosity(1:nLons, 1:nLats, 0:nAlts+1, iUp_) = CondResult

        VerticalVelocityStage1(1:nLons,1:nLats, 0:nAlts+1, iSpecies ) = &
             VerticalVelocity(1:nLons,1:nLats, 0:nAlts+1, iSpecies ,iBlock) + &
             Viscosity(1:nLons,1:nLats, 0:nAlts+1, iUp_  )

        DtCSLocal = Dt/2.0
        TmpVel = VerticalVelocityStage1(1:nLons, 1:nLats,-1:nAlts+2, iSpecies)
        TmpDiff = (4.0/3.0)*ViscCoefS(1:nLons, 1:nLats, 0:nAlts+1, iSpecies) + &
             Mass(iSpecies)*NDensityS(1:nLons, 1:nLats, 0:nAlts+1, iSpecies, iBlock) * &
             KappaEddyDiffusion(1:nLons, 1:nLats, 0:nAlts+1, iBlock)
        TmpMulFac = Mass(iSpecies)*NDensityS(1:nLons, 1:nLats, 0:nAlts+1, iSpecies, iBlock)

        call calc_conduction(iBlock, DtCSLocal, VertVelNeuBCS(iSpecies), &
             TmpVel, TmpDiff, TmpMulFac, CondResult)
        Viscosity(1:nLons, 1:nLats, 0:nAlts+1, iUp_) = CondResult

        VerticalVelocityH(1:nLons,1:nLats, 0:nAlts+1, iSpecies ) = &
             VerticalVelocityStage1(1:nLons,1:nLats, 0:nAlts+1, iSpecies ) + &
             Viscosity(1:nLons,1:nLats, 0:nAlts+1, iUp_  )

        ! Full Time Step Update
        DtCSLocal = Dt
        VelocityStage0 = VerticalVelocity(1:nLons, 1:nLats,-1:nAlts+2, iSpecies,iBlock)
        TmpDiff = (4.0/3.0)*ViscCoefS(1:nLons, 1:nLats, 0:nAlts+1, iSpecies) + &
             Mass(iSpecies)*NDensityS(1:nLons, 1:nLats, 0:nAlts+1, iSpecies, iBlock) * &
             KappaEddyDiffusion(1:nLons, 1:nLats, 0:nAlts+1, iBlock)
        TmpMulFac = Mass(iSpecies)*NDensityS(1:nLons, 1:nLats, 0:nAlts+1, iSpecies, iBlock)
        call calc_conduction(iBlock, DtCSLocal, VertVelNeuBCS(iSpecies), &
             VelocityStage0, TmpDiff, TmpMulFac, CondResult)
        Viscosity(1:nLons, 1:nLats, 0:nAlts+1, iUp_) = CondResult

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

     !BP (5/29/2020)                                                             
     !Compute the IR heating at every location in GITM.                             
     !Inputs:   SZA & Altitude                                                             
     !Outputs:  Q_IR (iLon, iLat, iAlt, iProc)                             
     QnirTOT(:,:,:,:) = 0.0
     if (UseIRHeating .and. isVenus) then
       call start_timing("IR Heating")
       do iAlt = 1, nAlts
         do iLon = 1, nLons
           do iLat = 1, nLats
             !Everywhere we don't want heating                                        
             !1.) Below 80 km                                                          
             !2.) Above 160 km                                                        
             !3.) SZA < 0                                                             
             !4.) SZA > 90                                                      
             if (Altitude_GB(iLon,iLat,iAlt,iBlock) < 80.0*1000.0 .or. &
                 Altitude_GB(iLon,iLat,iAlt,iBlock) > 160.0*1000.0 .or. &
                 sza(iLon,iLat,iBlock) < 0.0 .or. &
                 sza(iLon,iLat,iBlock) .ge. 89.0*pi/180.0) then
               QnirTOT(iLon,iLat,iAlt,iBlock) = 0.0
             else
               do iSza = 1, 11
                 if (sza(iLon, iLat, iBlock) .le. sza_table(iSza+1) .and. &
                     sza(iLon, iLat, iBlock) .ge. sza_table(iSza)) then
                   exit
                 endif
               enddo

               do jAlt = 1,40
                 if (Altitude_GB(iLon,iLat,iAlt,iBlock) .le. altitude_table(jAlt+1) .and. &
                     Altitude_GB(iLon,iLat,iAlt,iBlock) .ge. altitude_table(jAlt)) then
                   exit
                 endif
               enddo

               !Debugging stuff                                                      
               !if (iProc == 0) then                                                
                ! write(*,*) iLon, iLat, iAlt, nLons                                    
                 !write(*,*) "GITM SZA", sza(iLon,iLat, iBlock)                               
                 !write(*,*) "SZA Found in the table:", sza_table(i)                          
                 !write(*,*) "Altitude Found in table:", altitude_table(j)              
                 !write(*,*) "Longitude", Longitude(iLon, iBlock) * 180.0/pi             
                 !write(*,*) "Latitude", Latitude(iLat, iBlock) * 180.0/pi                
                 !write(*,*) "Altitude",  Altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0      
               !endif                                                  
               r_sza = (sza(iLon, iLat, iBlock) - sza_table(iSza)) / &
                       (sza_table(iSza+1) - sza_table(iSza))

               m = (altitude_table(jAlt) - altitude_table(jAlt+1))/&
                   (qIR_NLTE_table(jAlt,iSza) - qIR_NLTE_table(jAlt+1,iSza))

               x03 = (Altitude_GB(iLon,iLat,iAlt,iBlock) - &
                     (altitude_table(jAlt+1) - m*qIR_NLTE_table(jAlt+1,iSza))) / m
               m = (altitude_table(jAlt)- altitude_table(jAlt+1))/&
                   (qIR_NLTE_table(jAlt,iSza+1) - qIR_NLTE_table(jAlt+1,iSza+1))

               x12 = (Altitude_GB(iLon,iLat,iAlt,iBlock) - &
                     (altitude_table(jAlt+1) - m*qIR_NLTE_table(jAlt+1,iSza+1))) / m


               if (r_sza > 1.0) then
                 write(*,*) "r is too big...", iLon, iLat, iAlt
                 write(*,*) "GITM SZA:", sza(iLon, iLat, iBlock)
                 write(*,*) "Table SZA:", sza_table(iSza), sza_table(iSza+1)
               endif

               if (x03 < 0.0 .or. x12 < 0.0) then
                 write(*,*) "First interpolated value is negative...this is wrong."
                 write(*,*) "x03", x03
                 write(*,*) "x12", x12
               endif

               QnirTOT(iLon, iLat, iAlt, iBlock) = &
                 x03*(1 - r_sza) + x12*r_sza

               !Given in K/s. No need for rho or cp because that             
               !converts J/s -> K/s. Need TempUnit to normalize it to               
               !GITM units                                                        
               QnirTOT(iLon, iLat, iAlt, iBlock) = QnirTOT(iLon, iLat, iAlt, iBlock)  &
                                                   / TempUnit(iLon,iLat,iAlt)! &        
                                                   ! / Rho(iLon,iLat,iAlt, iBlock) &       
                                                   ! / cp(iLon,iLat,iAlt,iBlock)        
             endif

             !LTE IR Heating Section                                                    
             if (Altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0 < 100.0) then
               LST = mod(UTime/3600.0 + Longitude(iLon, iBlock)*180.0/(360*pi/RP_hours), &
                         RP_hours)

               !Convert to a 0 - 24 scale instead of the 0 - VenusHoursPerDay             
               LST = 24*LST/RP_hours

               !0 - 90 in 5 deg increments                            
               do iLatTable = 1,19
                 if (ABS(Latitude(iLat,iBlock)*180.0/pi) .ge. 5*(iLatTable - 1) .and. &
                     ABS(Latitude(iLat,iBlock)*180.0/pi) .le. 5*(iLatTable)) then
                   exit
                 endif
               enddo

               do iAltTable = 1,16
                 if (Altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0 .ge. 70 + 2*(iAltTable - 1) &
                     .and. &
                     Altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0 .le. 70 + 2*(iAltTable)) then
                   exit
                 endif
               enddo

               QnirTOT(iLon,iLat,iAlt,iBlock) = diurnalHeating(iLatTable,iAltTable)* &
                                                cos(2*pi/24*(LST - 12)) + &
                                                semidiurnalHeating(iLatTable,iAltTable)* &
                                                cos(4*pi/24*(LST - 12))

               if (QnirTOT(iLon,iLat,iAlt,iBlock) < 0.0) then
                 QnirTOT(iLon,iLat,iAlt,iBlock) = 0.0
               endif

               !Given in K/s. No need for rho or cp because that                              
               !converts J/s -> K/s. Need TempUnit to normalize it to              
               !GITM units                                                       
               QnirTOT(iLon, iLat, iAlt, iBlock) = QnirTOT(iLon, iLat, iAlt, iBlock)  &
                                                   / TempUnit(iLon,iLat,iAlt)! &      
                                                   ! / Rho(iLon,iLat,iAlt, iBlock) &       
                                                   ! / cp(iLon,iLat,iAlt,iBlock)           
             endif
           enddo
         enddo
       enddo
       call end_timing("IR Heating")
     endif
 
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
