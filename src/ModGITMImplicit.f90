!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModGITMImplicit

  !DESCRIPTION:
  ! This module implements a point implicit scheme for the implicit
  ! part of the right hand side Rimp = R - Rexp that contributes to the
  ! the implicitly treated variables Uimp, a subset of U=(Uexp, Uimp).
  ! Rimp should depend on the local cell values only: no spatial derivatives!
  !
  ! For a one stage scheme the variables are updated with the following
  ! first order scheme
  !
  ! Uexp^n+1 = Uexp^n + Dt*Rexp(U^n)
  ! Uimp^n+1 = Uimp^n + Dt*Rimp(U^n+1)
  !
  ! For the two stage scheme the following scheme is applied
  !
  ! Uexp^n+1/2 = Uexp^n + Dt/2 * Rexp(U^n)
  ! Uimp^n+1/2 = Uimp^n + Dt/2 * Rimp(U^n+1/2)
  !
  ! Uexp^n+1   = Uexp^n + Dt*Rexp(U^n+1/2)
  ! Uimp^n+1   = Uimp^n + Dt*beta*Rimp(U^n+1) + Dt*(1-beta)*Rimp(U^n)
  !
  ! where beta is in the range 0.5 and 1.0.
  ! The scheme is second order accurate in time for beta = 0.5.
  !
  ! For the general case Rimp is non-linear, and it is linearized as
  !
  ! Rimp(U^n+1/2) = Rimp(Uexp^n+1/2,Uimp^n) + dRimp/dUimp*(Uimp^n+1/2 - Uimp^n)
  ! Rimp(U^n+1)   = Rimp(Uexp^n+1,  Uimp^n) + dRimp/dUimp*(Uimp^n+1   - Uimp^n)
  !
  ! Note that the Jacobian dRimp/dUimp is evaluated at the partially advanced
  ! states (Uexp^n+1/2,Uimp^n) and (Uexp^n+1,Uimp^n) respectively.
  ! If Rimp is linear, the linearization is exact.
  !
  ! Substituting the linearization back into the one-stage and two-stage
  ! schemes yields a linear equation for the differences 
  ! (Uimp^n+1/2 - Uimp^n) and (Uimp^n+1 - Uimp^n), respectively.
  ! Since Rimp depends on the local cell values only, the linear equations
  ! can be solved point-wise for every cell.
  !
  ! The Jacobian can be given analytically by the subroutine passed to
  ! update_point_implicit, or it can be obtained by taking numerical
  ! derivatives of Rimp:
  !
  ! dRimp/dU^w = ((Rimp(Uexp^n+1,Uimp^n+eps^w) - Rimp(Uexp^n+1,Uimp^n))/eps^w
  !                
  ! where eps^w is a small perturbation in the w-th component of Uimp.
  !EOP

  implicit none

  save

  private ! except

  logical, public :: UsePointImplicit  ! Use point impl scheme?
  integer, public, allocatable :: &
       iVarPointImpl_I(:)                        ! Indexes of point impl. vars
  logical, public :: IsPointImplMatrixSet=.false.! Is dS/dU matrix analytic?
  logical, public :: IsPointImplPerturbed=.false.! Is the state perturbed?

  real, public, allocatable :: &
       DsDu_VVC(:,:,:,:,:), &     ! dS/dU derivative matrix
       EpsPointImpl_V(:)          ! absolute perturbation per variable
  real, public    :: EpsPointImpl ! relative perturbation

  public update_point_implicit    ! do update with point implicit scheme


  ! Local variables
  ! Number of point implicit variables
  integer :: nVarPointImpl   

contains
 
  !===========================================================================
  subroutine update_point_implicit(iLon,iLat,iAlt,iBlock,nSources,iSources,ChemicalHeating)

    use ModChemistry
    use ModGITM,  ONLY: iProc, dt,SpeciesDensity, SpeciesDensityOld
    use ModKind,    ONLY: nByteReal
    use ModPlanet, ONLY: nSpeciesAll
    use ModTime, ONLY: istep
    use ModInputs, ONLY: CFL,BetaPointImpl, IsASymmetric

    integer, intent(in) :: iBlock, ILon,iLat,iAlt
    real, intent(out) :: nSources(nSpeciesTotal),iSources(nions-1),chemicalheating

    real  :: IonSources(nIons),IonLosses(nIons)
    real  :: NeutralSources(nSpeciesTotal), NeutralLosses(nSpeciesTotal)
    real  :: ChemicalHeatingSub,Emission(nEmissions)

    integer :: nvar,iivar,ivar,ijvar,jvar, ispecies

    real :: DtCell, BetaStage, Norm, Epsilon
    real :: StateExpl_VC(nSpeciesAll),State_VGB(nSpeciesAll), StateOld_VCB(nSpeciesAll)
    real :: Source0_VC(nSpeciesAll), Source1_VC(nSpeciesAll),Source_VC(nSpeciesAll)
    real :: State0_C,DsDu_VVC(nSpeciesAll, nSpeciesAll),Matrix_II(nSpeciesAll,nSpeciesAll)
    real :: Rhs_I(nSpeciesAll)

    ! The default values for the state variables:
    ! Variables which are physically positive should be set to 1,
    ! variables that can be positive or negative should be set to 0:
    real, parameter :: DefaultState_V(nSpeciesAll) = 1

    character(len=*), parameter:: NameSub='update_point_implicit'
    character (len=100) :: test_string=''

    logical :: DoTest, DoTestMe,DoTestCell
    logical :: DoReplaceDensity = .true.
    !-------------------------------------------------------------------------

    call start_timing(NameSub)  
    nVar = nSpeciesAll
    ChemicalHeating = 0.0

    ! Initialization
    if(.not.allocated(iVarPointImpl_I))then

       ! Set default perturbation parameters
       allocate(EpsPointImpl_V(nVar))
       if(nByteReal == 8)then
          if(IsAsymmetric)then
             EpsPointImpl   = 1.e-6
          else 
             EpsPointImpl   = 1.e-9
          end if
          EpsPointImpl_V = 1.e-12
       else
          EpsPointImpl   = 1.e-3
          EpsPointImpl_V = 1.e-6
       end if

       ! This call should allocate and set the iVarPointImpl_I index array,
       ! set IsPointImplMatrixSet=.true. if the dS/dU matrix is analytic,
       ! it may also modify the EpsPointImpl and EpsPointImpl_V parameters.

       call init_pt_implicit(nvar)

       nVarPointImpl = size(iVarPointImpl_I)

       if(.not.allocated(iVarPointImpl_I)) call stop_GITM( &
            'calc_user_sources did not set iVarPointImpl_I')

    end if

    ! The beta parameter is always one in the first stage
    if(iStep == 1)then
       BetaStage = 1.0
    else
       BetaStage = BetaPointImpl
    end if

    ! Store explicit update
    StateExpl_VC(1:nSpeciesAll) = SpeciesDensity(iLon,iLat,iAlt,:,iBlock) ! After  the solver

    ! Put back old values into the implicit variables
    StateOld_VCB = SpeciesDensityOld(iLon,iLat,iAlt,:,iBlock) ! Before the solver
    State_VGB(1:nSpeciesAll) = StateOld_VCB

    ! Calculate unperturbed source for right hand side 
    ! and possibly also set analytic Jacobean matrix elements.
    ! Multi-ion may set its elements while the user uses numerical Jacobean.
    Source_VC = 0.0
    DsDu_VVC  = 0.0

    Neutrals = State_VGB(1:nSpeciesTotal)
    Ions(1:nIons-1) = State_VGB(nSpeciesTotal+1:nSpeciesAll)
    Ions(nions) = iDensityS(iLon,iLat,iAlt,nIons,iBlock)

!    call calc_reaction_rates(iLon,iLat,iAlt,iBlock)
    call calc_chemical_sources(iLon,iLat,iAlt,iBlock,IonSources,IonLosses,NeutralSources, &
         NeutralLosses,ChemicalHeatingSub,Emission)

    Source_VC(1:nSpeciesTotal) = NeutralSources - NeutralLosses
    Source_VC(nSpeciesTotal+1:nVar) = IonSources - IonLosses

    ! Calculate (part of) Jacobean numerically if necessary
    if(.not.IsPointImplMatrixSet)then
       ! Let the source subroutine know that the state is perturbed
       IsPointImplPerturbed = .true.

       ! Save unperturbed source
       Source0_VC = Source_VC(1:nVar)

       ! Perturb all point implicit variables one by one

       do iIVar = 1,nVar
          ivar = iivar

          ! Store unperturbed state
          State0_C = State_VGB(iVar)

          ! Get perturbation based value
          Norm = State0_C

          Epsilon = EpsPointImpl*Norm + EpsPointImpl_V(iVar)

          if(DefaultState_V(iVar) > 0.5 .and. .not. IsAsymmetric) &
               Epsilon = min(Epsilon, max(1e-30, 0.5*State0_C))

          ! Perturb the state
          State_VGB(iVar) = State0_C + Epsilon

          ! Calculate perturbed source
          Source_VC = 0.0

          !Must update the state in calc_chemical_sources!!!
          Neutrals = State_VGB(1:nSpeciesTotal)
          Ions(1:nIons-1) = State_VGB(nSpeciesTotal+1:nSpeciesAll)
          Ions(nions) = iDensityS(iLon,iLat,iAlt,nIons,iBlock)

          call calc_chemical_sources(iLon,iLat,iAlt,iBlock,IonSources,IonLosses,NeutralSources, &
               NeutralLosses,ChemicalHeatingSub,Emission) 


          Source_VC(1:nSpeciesTotal) = NeutralSources - NeutralLosses
          Source_VC(nSpeciesTotal+1:nVar) = IonSources - IonLosses

          if(IsAsymmetric)then

             ! Calculate dS/dU matrix elements
             do iJVar = 1,nVar
                jvar = ijvar
                DsDu_VVC(jVar,iVar) = DsDu_VVC(jVar,iVar) + &
                     (Source_VC(jVar) - Source0_VC(jVar))/Epsilon
             end do

          else

             ! Store perturbed source corresponding to +Epsilon perturbation
             Source1_VC = Source_VC(1:nVar)

             ! Perturb the state in opposite direction
             State_VGB(iVar) = State0_C - Epsilon

             ! Calculate perturbed source
             Source_VC = 0.0

             Neutrals = State_VGB(1:nSpeciesTotal)
             Ions(1:nIons-1) = State_VGB(nSpeciesTotal+1:nSpeciesAll)
             Ions(nions) = iDensityS(iLon,iLat,iAlt,nIons,iBlock)

             call calc_chemical_sources(iLon,iLat,iAlt,iBlock,IonSources,IonLosses,NeutralSources, &
                  NeutralLosses,ChemicalHeatingSub,Emission)

             Source_VC(1:nSpeciesTotal) = NeutralSources - NeutralLosses
             Source_VC(nSpeciesTotal+1:nVar) = IonSources - IonLosses

             ! Calculate dS/dU matrix elements with symmetric differencing
             do iJVar = 1,nVar
                jvar = ijvar
                DsDu_VVC(jVar,iVar) = DsDu_VVC(jVar,iVar) + &
                     0.5*(Source1_VC(jVar) - Source_VC(jVar)) &
                     /Epsilon
             end do

          end if

          !Restore unperturbed state
          State_VGB(iVar) = State0_C

       end do

       ! Restore unperturbed source
       Source_VC(1:nVar) = Source0_VC


       IsPointImplPerturbed = .false.

    end if

    ! Do the implicit update

    ! Do not update body cells

    DtCell = dt

    ! The right hand side is Uexpl - Uold + Sold
    do iIVar = 1, nVar
       ivar = iivar
       Rhs_I(iIVar) = StateExpl_VC(iVar) &
            - StateOld_VCB(iVar) &
            + DtCell * Source_VC(iVar)
    end do

    ! The matrix to be solved for is A = (I - beta*Dt*dS/dU)
    do iIVar = 1, nVar
       ivar = iivar
       do iJVar = 1, nVar
          jvar = ijvar
          Matrix_II(iIVar, iJVar) = - BetaStage*DtCell* &
               DsDu_VVC( iVar, jVar)
       end do
       ! Add unit matrix
       Matrix_II(iIVar,iIVar) = Matrix_II(iIVar,iIVar) + 1.0

    end do


    ! Solve the A.dU = RHS equation
    call linear_equation_solver(nVarPointImpl, Matrix_II, Rhs_I)

    ! Update: U^n+1 = U^n + dU
    do iIVar = 1, nVar
       ivar = iivar
       State_VGB(iVar) =&
            StateOld_VCB(iVar) + Rhs_I(iIVar)

    end do

    ! Fix negative species densities
    State_VGB = max(1e-5, State_VGB)

    nSources = State_VGB(1:nSpeciesTotal) - StateExpl_VC(1:nSpeciesTotal)
    iSources = State_VGB(nSpeciesTotal+1:nVar) - StateExpl_VC(nSpeciesTotal+1:nvar)

!  if (ialt .eq. 42) then
!     write(*,*) "impl:"
!     do ivar = 1, nspeciestotal
!        write(*,*) ialt,"end chem: ",ivar,state_VGB(ivar),stateexpl_VC(ivar),nsources(ivar)
!     enddo
!     do ivar = 1, nions - 1
!        write(*,*) ialt,"end chem: ",ivar,state_vgb(ivar+nspeciestotal),&
!             stateexpl_vc(ivar+nspeciestotal),isources(ivar)
!     enddo
!
!     stop
!     endif

    call end_timing(NameSub)

  end subroutine update_point_implicit

  !============================================================================

  subroutine linear_equation_solver(nVar, Matrix_VV, Rhs_V)

    integer, intent(in) :: nVar
    real, intent(inout) :: Matrix_VV(nVar, nVar)
    real, intent(inout) :: Rhs_V(nVar)

    ! This routine solves the system of Nvar linear equations:
    ! 
    !               Matrix_VV*dUCell = Rhs_V.
    ! 
    ! The result is returned in Rhs_V, the matrix is overwritten
    ! with the LU decomposition.
    !
    ! The routine performs a lower-upper (LU) decomposition of the 
    ! square matrix Matrix_VV of rank Nvar and then uses forward and
    ! backward substitution to obtain the solution vector dUCell.
    ! Crout's method with partial implicit pivoting is used to perform
    ! the decompostion.

    integer, parameter :: MAXVAR = 100

    integer :: IL, II, ILMAX, JL, KL, LL, INDX(MAXVAR)
    real    :: SCALING(MAXVAR), LHSMAX, LHSTEMP, TOTALSUM
    real, parameter :: TINY=1.0E-20

    !--------------------------------------------------------------------------
    if(nVar > MAXVAR) call stop_GITM(&
         'ERROR in ModPointImplicit linear solver: MaxVar is too small')

    !\
    ! Loop through each row to get implicit scaling
    ! information.
    !/
    DO IL=1,nVar
       LHSMAX=0.00
       DO JL=1,nVar
          IF (ABS(Matrix_VV(IL,JL)).GT.LHSMAX) LHSMAX=ABS(Matrix_VV(IL,JL))
       END DO
       SCALING(IL)=1.00/LHSMAX
    END DO

    !\
    ! Peform the LU decompostion using Crout's method.
    !/
    DO JL=1,nVar
       DO IL=1,JL-1
          TOTALSUM=Matrix_VV(IL,JL)
          DO KL=1,IL-1
             TOTALSUM=TOTALSUM-Matrix_VV(IL,KL)*Matrix_VV(KL,JL)
          END DO
          Matrix_VV(IL,JL)=TOTALSUM
       END DO
       LHSMAX=0.00
       DO IL=JL,nVar
          TOTALSUM=Matrix_VV(IL,JL)
          DO KL=1,JL-1
             TOTALSUM=TOTALSUM-Matrix_VV(IL,KL)*Matrix_VV(KL,JL)
          END DO
          Matrix_VV(IL,JL)=TOTALSUM
          LHSTEMP=SCALING(IL)*ABS(TOTALSUM)
          IF (LHSTEMP.GE.LHSMAX) THEN
             ILMAX=IL
             LHSMAX=LHSTEMP
          END IF
       END DO
       IF (JL.NE.ILMAX) THEN
          DO KL=1,nVar
             LHSTEMP=Matrix_VV(ILMAX,KL)
             Matrix_VV(ILMAX,KL)=Matrix_VV(JL,KL)
             Matrix_VV(JL,KL)=LHSTEMP
          END DO
          SCALING(ILMAX)=SCALING(JL)
       END IF
       INDX(JL)=ILMAX
       IF (abs(Matrix_VV(JL,JL)).EQ.0.00) Matrix_VV(JL,JL)=TINY
       IF (JL.NE.nVar) THEN
          LHSTEMP=1.00/Matrix_VV(JL,JL)
          DO IL=JL+1,nVar
             Matrix_VV(IL,JL)=Matrix_VV(IL,JL)*LHSTEMP
          END DO
       END IF
    END DO

    !\
    ! Peform the forward and back substitution to obtain
    ! the solution vector.
    !/
    II=0
    DO IL=1,nVar
       LL=INDX(IL)
       TOTALSUM=Rhs_V(LL)
       Rhs_V(LL)=Rhs_V(IL)
       IF (II.NE.0) THEN
          DO JL=II,IL-1
             TOTALSUM=TOTALSUM-Matrix_VV(IL,JL)*Rhs_V(JL)
          END DO
       ELSE IF (TOTALSUM.NE.0.00) THEN
          II=IL
       END IF
       Rhs_V(IL)=TOTALSUM
    END DO
    DO IL=nVar,1,-1
       TOTALSUM=Rhs_V(IL)
       DO JL=IL+1,nVar
          TOTALSUM=TOTALSUM-Matrix_VV(IL,JL)*Rhs_V(JL)
       END DO
       Rhs_V(IL)=TOTALSUM/Matrix_VV(IL,IL)
    END DO

    
  end subroutine linear_equation_solver
  
  subroutine init_pt_implicit(nvar)

    integer, intent(in) :: nvar

    integer :: iVar
    
    if(allocated(iVarPointImpl_I)) deallocate(iVarPointImpl_I)
    allocate(iVarPointImpl_I(nVar))
    
    do iVar = 1, nVar
       iVarPointImpl_I(iVar) = iVar
    enddo

  end subroutine init_pt_implicit
    

end module ModGITMImplicit
