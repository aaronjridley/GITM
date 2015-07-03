!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!\
! ------------------------------------------------------------
! advance
! ------------------------------------------------------------
!/

subroutine advance_vertical_1d

  use ModVertical
  use ModGITM, ONLY : Dt, iCommGITM, iProc
  use ModInputs, only: UseBarriers, iDebugLevel
  implicit none

  !-----------------------------------------------------------

  integer :: iError

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 6) write(*,*) "=======> vertical bcs 1", iproc

  ! Fill in ghost cells
  call set_vertical_bcs(LogRho,LogNS,Vel_GD,Temp,LogINS,IVel,VertVel)

  ! Copy input state into New state
  NewLogNS  = LogNS
  NewLogINS = LogINS
  NewLogRho = LogRho
  NewVel_GD = Vel_GD
  NewTemp   = Temp
  NewVertVel = VertVel

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 7) write(*,*) "========> stage 1", iproc

  ! Do the half step: U^n+1/2 = U^n + (Dt/2) * R(U^n)
!  Dt = Dt/2
!
!  write(*,*) "vv, before av1s 1: ", VertVel(3, :),LogNS(3,1)
!
!  call advance_vertical_1stage(&
!       LogRho, LogNS, Vel_GD, Temp, NewLogRho, NewLogNS, NewVel_GD, NewTemp, &
!       LogINS, NewLogINS, IVel, VertVel, NewVertVel)
!
!  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
!  if (iDebugLevel > 7) write(*,*) "========> vertical bcs 2", iproc
!
!  ! Fill in ghost cells for U^n+1/2 state
!  call set_vertical_bcs(NewLogRho, NewLogNS, NewVel_GD, NewTemp,NewLogINS, &
!       IVel, NewVertVel)
!
!  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
!  if (iDebugLevel > 7) write(*,*) "========> stage 2", iproc
!
!  write(*,*) "vv, before av1s 2: ", NewVertVel(3, :),NewLogNS(3,1)
!
!  ! Do full step U^n+1 = U^n + Dt * R(U^n+1/2)
!  Dt = 2*Dt
  call advance_vertical_1stage(&
       NewLogRho, NewLogNS, NewVel_GD, NewTemp, LogRho, LogNS, Vel_GD, Temp, &
       NewLogINS, LogINS, IVel, NewVertVel, VertVel)
  
!  write(*,*) "vv, after av1s 2: ", VertVel(3, :),LogNS(3,1)

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 7) write(*,*) "========> vertical bcs 3", iproc

  ! Fill in ghost cells for updated U^n+1 state
  call set_vertical_bcs(LogRho, LogNS, Vel_GD, Temp, LogINS, IVel, VertVel)

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 7) &
       write(*,*) "========> Done with advance_vertical_1d", iproc

end subroutine advance_vertical_1d

!=============================================================================
subroutine advance_vertical_1stage( &
     LogRho, LogNS, Vel_GD, Temp, NewLogRho, NewLogNS, NewVel_GD, NewTemp, &
     LogINS, NewLogINS, IVel, VertVel, NewVertVel)

  ! With fluxes and sources based on LogRho..Temp, update NewLogRho..NewTemp

  use ModGITM, only: &
       Dt, iEast_, iNorth_, iUp_
  use ModPlanet
  use ModSizeGitm
  use ModVertical, only : &
       Heating, EddyCoef_1D, ViscCoef_1d,Centrifugal, Coriolis, &
       MeanMajorMass_1d, Gamma_1d, InvRadialDistance_C, &
       Gravity_G, Altitude_G,Cv_1D
  use ModTime
  use ModInputs
  use ModConstants
  use ModSources, only : EddyCondAdia
  use ModUserGITM, only: userdata1d

  implicit none

  real, intent(in) :: LogRho(-1:nAlts+2)
  real, intent(in) :: LogNS(-1:nAlts+2,nSpecies)
  real, intent(in) :: LogINS(-1:nAlts+2,nIonsAdvect)
  real, intent(in) :: Vel_GD(-1:nAlts+2,3)
  real, intent(in) :: IVel(-1:nAlts+2,3)
  real, intent(in) :: Temp(-1:nAlts+2)
  real, intent(in) :: VertVel(-1:nAlts+2,nSpecies)

  real, intent(inout) :: NewLogRho(-1:nAlts+2)
  real, intent(inout) :: NewLogNS(-1:nAlts+2,nSpecies)
  real, intent(inout) :: NewLogINS(-1:nAlts+2,nIonsAdvect)
  real, intent(inout) :: NewVel_GD(-1:nAlts+2,3)
  real :: NewVel2_G(-1:nAlts+2)
  real, intent(inout) :: NewTemp(-1:nAlts+2)
  real, intent(out) :: NewVertVel(-1:nAlts+2,nSpecies)
  real :: NS(-1:nAlts+2,nSpecies), Pressure1D(-1:nAlts+2)
  real :: Rho(-1:nAlts+2)

  real :: LogNum(-1:nAlts+2)

  real, dimension(1:nAlts)    :: GradLogRho, DivVel, GradTemp, GradTempKoM, &
       DiffLogRho, DiffTemp, GradTmp, DiffTmp, DiffLogNum, GradLogNum
  real, dimension(1:nAlts,3) :: GradVel_CD, DiffVel_CD

  real, dimension(1:nAlts,nSpecies)    :: GradLogNS, DiffLogNS, &
       GradVertVel, DiffVertVel, DivVertVel
  real, dimension(1:nAlts,nIonsAdvect) :: GradLogINS, DiffLogINS
  real :: NewSumRho, NewLogSumRho, rat, ed

  integer :: iAlt, iSpecies, jSpecies, iDim

  real, dimension(-1:nAlts+2)    :: NT
  real, dimension(-1:nAlts+2)    :: Press, LogPress
  real, dimension(1:nAlts)    :: DiffLogPress, GradLogPress
  real, dimension(1:nAlts,nSpecies)    :: EddyDiffusionVel

  real :: nVel(1:nAlts,1:nSpecies)
  integer :: nFilter, iFilter
  real :: LowFilter

!! WAVEDRAG Heating  Hickey et al [2000]
  real, dimension(1:nAlts)    :: StressHeating

!\
! Parameters Used for the Sponge
! This Sponge is useful to dampen out spurious modes
! oscillating between the bottom and top of the model.

  integer :: nAltsSponge = 12
  real :: kSP, NuSP, AmpSP



  !--------------------------------------------------------------------------

  NS = exp(LogNS)
  Rho = exp(LogRho)
  LogNum = alog(sum(NS,dim=2))
  nFilter = 10
  

     NT(-1:nAlts+2) = exp(LogNum(-1:nAlts+2))
  do iAlt = -1, nAlts + 2
    Press(iAlt) = NT(iAlt)*Boltzmanns_Constant*Temp(iAlt)
    LogPress(iAlt) = alog(Press(iAlt))
  enddo

  call calc_rusanov_alts(LogPress ,GradLogPress,  DiffLogPress)
  call calc_rusanov_alts(LogRho ,GradLogRho,  DiffLogRho)
  call calc_rusanov_alts(LogNum ,GradLogNum,  DiffLogNum)
  call calc_rusanov_alts(Temp   ,GradTemp,    DiffTemp)
  do iDim = 1, 3
     call calc_rusanov_alts(Vel_GD(:,iDim), &
          GradVel_CD(:,iDim),DiffVel_CD(:,iDim))
  enddo

  ! Add geometrical correction to gradient and obtain divergence
  DivVel = GradVel_CD(:,iUp_) + 2*Vel_GD(1:nAlts,iUp_)*InvRadialDistance_C

  do iSpecies=1,nSpecies

     call calc_rusanov_alts(LogNS(:,iSpecies),GradTmp, DiffTmp)
     GradLogNS(:,iSpecies) = GradTmp
     DiffLogNS(:,iSpecies) = DiffTmp

     call calc_rusanov_alts(VertVel(:,iSpecies),GradTmp, DiffTmp)
     GradVertVel(:,iSpecies) = GradTmp
     DiffVertVel(:,iSpecies) = DiffTmp
     DivVertVel(:,iSpecies) = GradVertVel(:,iSpecies) + &
          2*VertVel(1:nAlts,iSpecies)*InvRadialDistance_C

  enddo

  do iSpecies=1,nIonsAdvect
     call calc_rusanov_alts(LogINS(:,iSpecies), GradTmp, DiffTmp)
     GradLogINS(:,iSpecies) = GradTmp
     DiffLogINS(:,iSpecies) = DiffTmp
  enddo


  AmpSP = (1.0/(10.0*Dt))
  kSP = nAltsSponge + 1

  do iAlt = 1,nAlts

     NewLogRho(iAlt) = NewLogRho(iAlt) - Dt * &
          (DivVel(iAlt) + Vel_GD(iAlt,iUp_) * GradLogRho(iAlt) ) &
          + Dt * DiffLogRho(iAlt)

     do iSpecies=1,nSpecies
        NewLogNS(iAlt,iSpecies) = LogNS(iAlt,iSpecies) - Dt * &
             (DivVertVel(iAlt,iSpecies) + &
             VertVel(iAlt,iSpecies) * GradLogNS(iAlt,iSpecies) ) + &
              Dt * DiffLogNS(iAlt,iSpecies)

     enddo

     do iSpecies=1,nIonsAdvect
        NewLogINS(iAlt,iSpecies) = NewLogINS(iAlt,iSpecies) - Dt * &
             (IVel(iAlt,iUp_) * GradLogINS(iAlt,iSpecies) ) &
             + Dt * DiffLogINS(iAlt,iSpecies)
     enddo

!     ! dVr/dt = -[ (V grad V)_r + grad T + T grad ln Rho - g ]
!     ! and V grad V contains the centripetal acceleration 
!     ! (Vphi**2+Vtheta**2)/R
!     NewVel_GD(iAlt,iUp_) = NewVel_GD(iAlt,iUp_) - Dt * &
!          (Vel_GD(iAlt,iUp_)*GradVel_CD(iAlt,iUp_) &
!          - (Vel_GD(iAlt,iNorth_)**2 + Vel_GD(iAlt,iEast_)**2) &
!          * InvRadialDistance_C(iAlt) &
!          - Gravity_G(iAlt)) &
!          + Dt * DiffVel_CD(iAlt,iUp_)

     NewVel_GD(iAlt,iUp_) = 0.0


     if (iAlt >= (nAlts - nAltsSponge)) then
 
      NuSP = AmpSP*(1.0 - cos( pi*(kSP - (nAlts - iAlt))/kSP) )

     else

       NuSP = 0.0

     endif

     do iSpecies=1,nSpecies

        ! Version of vertical velocity with grad(p) and g here :

        NewVertVel(iAlt, iSpecies) = VertVel(iAlt, iSpecies) - Dt * &
             (VertVel(iAlt,iSpecies)*GradVertVel(iAlt,iSpecies) &
             - (Vel_GD(iAlt,iNorth_)**2 + Vel_GD(iAlt,iEast_)**2) &
             * InvRadialDistance_C(iAlt) + &
             Temp(iAlt)*GradLogNS(iAlt,iSpecies) * Boltzmanns_Constant / &
             Mass(iSpecies) + &
             GradTemp(iAlt) * Boltzmanns_Constant / Mass(iSpecies) &
             - Gravity_G(iAlt)) &
             + Dt * DiffVertVel(iAlt,iSpecies) 

        
        if (UseCoriolis) then
           NewVertVel(iAlt,ispecies) = NewVertVel(iAlt,ispecies) + Dt * ( &
                Centrifugal / InvRadialDistance_C(iAlt) + &
                Coriolis * Vel_GD(iAlt,iEast_))
        endif

     enddo
    
  enddo

  do iAlt = 1, nAlts

     do iSpecies=1,nSpecies

        NewVertVel(iAlt, iSpecies) = max(-MaximumVerticalVelocity, NewVertVel(iAlt, iSpecies))
        NewVertVel(iAlt, iSpecies) = min( MaximumVerticalVelocity, NewVertVel(iAlt, iSpecies))

        NewVel_GD(iAlt,iUp_) = NewVel_GD(iAlt,iUp_) + &
             NewVertVel(iAlt, iSpecies) * &
             (Mass(iSpecies) * NS(iAlt,iSpecies) / Rho(iAlt))

     enddo

  enddo

  StressHeating = 0.0

  if (UseStressHeating) then

    do iAlt = 1, nAlts 

      StressHeating(iAlt) = ViscCoef_1d(iAlt)* &
       (  (  (Gamma_1d(iAlt) - 1.0)/ ( NT(iAlt)*Boltzmanns_Constant) ) * &
           (  &
              (4.0/3.0)*GradVel_CD(iAlt,iUp_)**2 +    &
                        GradVel_CD(iAlt,iNorth_)**2 + &
                        GradVel_CD(iAlt,iEast_)**2    &
           )  )

    enddo

  endif
  userdata1d = 0.0
  do iAlt = 1, nAlts

     ! dVphi/dt = - (V grad V)_phi
     NewVel_GD(iAlt,iEast_) = NewVel_GD(iAlt,iEast_) - Dt * &
          Vel_GD(iAlt,iUp_)*GradVel_CD(iAlt,iEast_) &
          + Dt * DiffVel_CD(iAlt,iEast_)

     ! dVtheta/dt = - (V grad V)_theta
     NewVel_GD(iAlt,iNorth_) = NewVel_GD(iAlt,iNorth_) - Dt * &
          Vel_GD(iAlt,iUp_)*GradVel_CD(iAlt,iNorth_) &
          + Dt * DiffVel_CD(iAlt,iNorth_)

     ! dT/dt = -(V.grad T + (gamma - 1) T div V +  &
     !        (gamma - 1) * g  * grad (KeH^2  * rho) /rho 

!     if (UseTurbulentCond) then
!        NewTemp(iAlt)   = NewTemp(iAlt) - Dt * &
!             (Vel_GD(iAlt,iUp_)*GradTemp(iAlt) + &
!             (Gamma_1d(iAlt) - 1.0) * Temp(iAlt)*DivVel(iAlt))&
!             + Dt * DiffTemp(iAlt)
!     else

        NewTemp(iAlt)   = NewTemp(iAlt) - Dt * &
             (Vel_GD(iAlt,iUp_)*GradTemp(iAlt) + &
             (Gamma_1d(iAlt) - 1.0) * ( &
             Temp(iAlt)*DivVel(iAlt))) &
             + Dt * DiffTemp(iAlt) & 
             + Dt * StressHeating(iAlt) 


userdata1D(1,1,ialt,2) = -Dt * (Vel_GD(iAlt,iUp_)*GradTemp(iAlt))
userdata1D(1,1,ialt,3) = -Dt * (Gamma_1d(iAlt) - 1.0) * ( &
             Temp(iAlt)*DivVel(iAlt))
userdata1D(1,1,ialt,4) =  Dt * DiffTemp(iAlt) 
userdata1D(1,1,ialt,5) =  Dt * StressHeating(iAlt) 

!        NewTemp(iAlt)   = NewTemp(iAlt) + &
!             Dt*Vel_GD(iAlt,iUp_)* &
!             ( (Gamma_1d(iAlt) - 1.0)/Gamma_1d(iAlt) )*&
!               Temp(iAlt)*GradLogPress(iAlt)

  end do

  do iAlt = 1, nAlts
     NewSumRho    = sum( Mass(1:nSpecies)*exp(NewLogNS(iAlt,1:nSpecies)) )
     NewLogRho(iAlt) = alog(NewSumRho)
  enddo

end subroutine advance_vertical_1stage

!\
! ------------------------------------------------------------
! calc_rusanov
! ------------------------------------------------------------
!/

subroutine calc_rusanov_alts(Var, GradVar, DiffVar)

  use ModSizeGitm
  use ModVertical, only : dAlt_C, cMax
  implicit none

  real, intent(in) :: Var(-1:nAlts+2)
  real, intent(out):: GradVar(1:nAlts), DiffVar(1:nAlts)

  real, dimension(1:nAlts+1) :: VarLeft, VarRight, DiffFlux
  !------------------------------------------------------------

  call calc_facevalues_alts(Var, VarLeft, VarRight)

  ! Gradient based on averaged Left/Right values
  GradVar = 0.5 * &
       (VarLeft(2:nAlts+1)+VarRight(2:nAlts+1) - &
       VarLeft(1:nAlts)-VarRight(1:nAlts))/dAlt_C

  ! Rusanov/Lax-Friedrichs diffusive term
  DiffFlux = 0.5 * max(cMax(0:nAlts),cMax(1:nAlts+1)) * (VarRight - VarLeft)

  DiffVar = (DiffFlux(2:nAlts+1) - DiffFlux(1:nAlts))/dAlt_C

end subroutine calc_rusanov_alts

!\
! ------------------------------------------------------------
! calc_facevalues_alts
! ------------------------------------------------------------
!/

subroutine calc_facevalues_alts(Var, VarLeft, VarRight)

  use ModVertical, only: dAlt_F, InvDAlt_F
  use ModSizeGITM, only: nAlts
  use ModLimiterGitm

  implicit none
  
  real, intent(in) :: Var(-1:nAlts+2)
  real, intent(out):: VarLeft(1:nAlts+1), VarRight(1:nAlts+1)

  real :: dVarUp, dVarDown, dVarLimited(0:nAlts+1)

  real, parameter :: Factor1=0.6250000 ! 15/24
  real, parameter :: Factor2=0.0416667 !  1/24
  real :: h

  integer :: i

  do i=1,nAlts

     ! 4th order scheme for calculating face values

     h  = InvDAlt_F(i+1)*2.0
     dVarUp   = h*(Factor1*(Var(i+1)-Var(i)  ) - Factor2*(Var(i+2)-Var(i-1)))
     h  = InvDAlt_F(i)*2.0
     dVarDown = h*(Factor1*(Var(i)  -Var(i-1)) - Factor2*(Var(i+1)-Var(i-2)))

!     ! This is Gabor's scheme
!     dVarUp            = (Var(i+1) - Var(i))   * InvDAlt_F(i+1)
!     dVarDown          = (Var(i)   - Var(i-1)) * InvDAlt_F(i)

     dVarLimited(i) = Limiter_mc(dVarUp, dVarDown)

!     write(*,*) dVarUp, dVarDown, dVarLimited(i)

  end do

  i = 0
  dVarUp            = (Var(i+1) - Var(i))   * InvDAlt_F(i+1)
  dVarDown          = (Var(i)   - Var(i-1)) * InvDAlt_F(i)
  dVarLimited(i) = Limiter_mc(dVarUp, dVarDown)

  i = nAlts+1
  dVarUp            = (Var(i+1) - Var(i))   * InvDAlt_F(i+1)
  dVarDown          = (Var(i)   - Var(i-1)) * InvDAlt_F(i)
  dVarLimited(i) = Limiter_mc(dVarUp, dVarDown)

  do i=1,nAlts+1
     VarLeft(i)  = Var(i-1) + 0.5*dVarLimited(i-1) * dAlt_F(i)
     VarRight(i) = Var(i)   - 0.5*dVarLimited(i)   * dAlt_F(i)
  end do

end subroutine calc_facevalues_alts


