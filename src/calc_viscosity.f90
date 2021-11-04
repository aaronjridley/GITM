! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

subroutine calc_viscosity(iBlock)

  use ModInputs
  use ModSources
  use ModGITM
  implicit none

  integer, intent(in) :: iBlock
  integer :: iErr, iDir, iSpecies

  ! NeuBCS used for Conduction and Bulk Winds
  logical ::NeuBCS                      ! Set to true if you want Neumann BCs
  ! Species-Specific Boundary Condition Settings
  logical ::VertVelNeuBCS(1:nSpecies)   ! Set to true if you want Neumann BCs

  ! Sub-timestep used in the multi-step conduction
  real :: DtCSLocal

  real :: TmpVel(nLons, nLats, -1:nAlts+2)
  real :: TmpDiff(nLons, nLats, 0:nAlts+1)
  real :: TmpMulFac(nLons, nLats, 0:nAlts+1)
  real :: Rho110(nLons, nLats, 0:nAlts+1)
  real :: CondResult(nLons, nLats, 0:nAlts+1)
  
  ! Bulk Viscosity Variables
  real :: VelocityStage0(1:nLons,1:nLats,-1:nAlts+2)
  real :: VelocityStage1(1:nLons,1:nLats,-1:nAlts+2,1:3)
  real :: VelocityH(1:nLons,1:nLats,-1:nAlts+2,1:3)
  real :: VelocityF(1:nLons,1:nLats,-1:nAlts+2,1:3)

  ! Vertical Viscosity Variables
  real :: VerticalVelocityStage0(1:nLons,1:nLats,-1:nAlts+2,1:nSpecies)
  real :: VerticalVelocityStage1(1:nLons,1:nLats,-1:nAlts+2,1:nSpecies)
  real :: VerticalVelocityH(1:nLons,1:nLats,-1:nAlts+2,1:nSpecies)
  real :: VerticalVelocityF(1:nLons,1:nLats,-1:nAlts+2,1:nSpecies)

  if (iDebugLevel > 4) write(*,*) "=====> Calc Viscosity", iproc
  if (UseBarriers) call MPI_BARRIER(iCommGITM, iErr)

  ! Note, this assumes that all species need Neumann Upper boundary
  ! conditions.  If you specify escape rates in set_vertical_bcs.f90
  ! such as Jeans Escape, then set VertVelNeuBCs(iH_ or iH2_) = .false.
  ! This will allow you to specify escape rates at the top.
  VertVelNeuBCS(1:nSpecies) = .true. 
  
  call calc_viscosity_coef(iBlock)

  ! 07/13/2017:  JMB, Below is the 2nd Order
  ! Implicit method for Horizontal Winds.
  ! We assume that these winds (North, East)
  ! have a zero gradient at the top, so NeuBCS = .true.

  Rho110 = Rho(1:nLons, 1:nLats,0:nAlts+1, iBlock)
     
  NeuBCS = .true.
  do iDir = iEast_, iNorth_

     if (iDebugLevel > 4) write(*,*) "=====> Viscosity stage 1", iproc, iDir, 1
     if (UseBarriers) call MPI_BARRIER(iCommGITM,iErr)

     DtCSLocal = Dt/2.0
     VelocityStage0 = Velocity(1:nLons, 1:nLats,-1:nAlts+2, iDir,iBlock)
     TmpDiff = ViscCoef(1:nLons, 1:nLats, 0:nAlts+1) + &
          Rho110 * KappaEddyDiffusion(1:nLons, 1:nLats, 0:nAlts+1, iBlock)

     call calc_conduction_1d(iBlock, DtCSLocal, NeuBCS, &
          VelocityStage0, TmpDiff, Rho110, CondResult)

     Viscosity(1:nLons, 1:nLats, 0:nAlts+1, iDir) = CondResult

     VelocityStage1(1:nLons,1:nLats, 0:nAlts+1, iDir ) = &
          Velocity(1:nLons,1:nLats, 0:nAlts+1, iDir ,iBlock) + &
          Viscosity(1:nLons,1:nLats, 0:nAlts+1, iDir  )

     if (iDebugLevel > 4) write(*,*) "=====> Viscosity stage 2", iproc, iDir, 2
     if (UseBarriers) call MPI_BARRIER(iCommGITM,iErr)

     DtCSLocal = Dt/2.0
     TmpVel = VelocityStage1(1:nLons, 1:nLats,-1:nAlts+2, iDir)
     TmpDiff = ViscCoef(1:nLons, 1:nLats, 0:nAlts+1) + &
          Rho110 * KappaEddyDiffusion(1:nLons, 1:nLats, 0:nAlts+1, iBlock)

     call calc_conduction_1d(iBlock, DtCSLocal, NeuBCS, &
          TmpVel, TmpDiff, Rho110, CondResult)

     Viscosity(1:nLons, 1:nLats, 0:nAlts+1,iDir) = CondResult

     VelocityH(1:nLons,1:nLats, 0:nAlts+1, iDir ) = &
          VelocityStage1(1:nLons,1:nLats, 0:nAlts+1, iDir ) + &
          Viscosity(1:nLons,1:nLats, 0:nAlts+1, iDir  )

     if (iDebugLevel > 4) write(*,*) "=====> Viscosity stage 3", iproc, iDir, 3
     if (UseBarriers) call MPI_BARRIER(iCommGITM,iErr)

     DtCSLocal = Dt
     VelocityStage0 = Velocity(1:nLons, 1:nLats,-1:nAlts+2, iDir,iBlock)
     TmpDiff = ViscCoef(1:nLons, 1:nLats, 0:nAlts+1) + &
          Rho110 * KappaEddyDiffusion(1:nLons, 1:nLats, 0:nAlts+1, iBlock)
     call calc_conduction_1d(iBlock, DtCSLocal, NeuBCS, &
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
  if (UseBarriers) call MPI_BARRIER(iCommGITM,iErr)

  ! JMB:  07/13/2017
  ! Add 2nd Order Viscosity for Each Species
  Viscosity(:,:,:,iUp_) = 0.0

  do iSpecies = 1, nSpecies

     DtCSLocal = Dt/2.0
     VelocityStage0 = VerticalVelocity(1:nLons, 1:nLats,-1:nAlts+2, iSpecies,iBlock)
     TmpDiff = (4.0/3.0) * ViscCoefS(1:nLons, 1:nLats, 0:nAlts+1, iSpecies) + &
          Mass(iSpecies) * NDensityS(1:nLons, 1:nLats, 0:nAlts+1, iSpecies, iBlock) * &
          KappaEddyDiffusion(1:nLons, 1:nLats, 0:nAlts+1, iBlock)
     TmpMulFac = Mass(iSpecies) * NDensityS(1:nLons, 1:nLats, 0:nAlts+1, iSpecies, iBlock)

     call calc_conduction_1d(iBlock, DtCSLocal, VertVelNeuBCS(iSpecies), &
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

     call calc_conduction_1d(iBlock, DtCSLocal, VertVelNeuBCS(iSpecies), &
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
     call calc_conduction_1d(iBlock, DtCSLocal, VertVelNeuBCS(iSpecies), &
          VelocityStage0, TmpDiff, TmpMulFac, CondResult)
     Viscosity(1:nLons, 1:nLats, 0:nAlts+1, iUp_) = CondResult

     VerticalVelocityF(1:nLons,1:nLats, 0:nAlts+1, iSpecies) = &
          VerticalVelocity(1:nLons,1:nLats, 0:nAlts+1, iSpecies,iBlock) + &
          Viscosity(1:nLons,1:nLats, 0:nAlts+1, iUp_  )

     !! Simply update the Vertical Winds here (rather than in add_sources)
     !VerticalVelocity(1:nLons,1:nLats, 0:nAlts+1, iSpecies, iBlock) = &
     !     2.0*VerticalVelocityH(1:nLons,1:nLats, 0:nAlts+1, iSpecies) - &
     !     VerticalVelocityF(1:nLons,1:nLats, 0:nAlts+1, iSpecies) 

     ! This is the final 2nd ORder Implicit Viscosity for 
     ! the horizontal Bulk Winds
     VerticalViscosityS(1:nLons, 1:nLats, 0:nAlts+1, iSpecies ) = &
          ( 2.0*VerticalVelocityH(1:nLons, 1:nLats, 0:nAlts+1, iSpecies) - &
          VerticalVelocityF(1:nLons, 1:nLats, 0:nAlts+1, iSpecies)) - &
          VerticalVelocity(1:nLons, 1:nLats, 0:nAlts+1, iSpecies, iBlock) 

  enddo ! iSpecies

end subroutine calc_viscosity
