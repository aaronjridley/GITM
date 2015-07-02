!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModHypre

  use ModPotentialField
  use ModKind, ONLY: Int8_

  implicit none

  private
  public:: hypre_initialize
  public:: hypre_solver
  public:: hypre_preconditioner
  public:: hypre_read_param

  ! These are defined in HYPREf.h (Fortran header file)
  integer, parameter:: HYPRE_STRUCT = 1111, HYPRE_PARCSR = 5555

  ! This is defined in HYPRE_sstruct_mv.h (a C header file)
  integer, parameter:: HYPRE_SSTRUCT_VARIABLE_CELL = 0

  ! We can use a single part with many boxes or many parts with local cell indexes
  logical:: UseSinglePart = .true.

  integer, parameter:: nDim = 3
  integer, parameter:: nStencil = 2*nDim+1

  ! Hypre uses C type 8-byte integers for pointers
  integer(Int8_):: i8Grid, i8Graph, i8Stencil, i8Precond, i8Solver
  integer(Int8_):: i8A, i8B, i8X, i8ParA, i8ParB, i8ParX

  integer:: nPart = 1 ! grid consists of one or more parts
  integer:: iPart = 0 ! index of the part
  integer:: jPart = 0 ! index of the neighbor part
  integer:: nVar  = 1 ! there is only one variable per grid cell
  integer:: iVar  = 0 ! index of the variable
  integer:: iVarType_I(1) = (/ HYPRE_SSTRUCT_VARIABLE_CELL /)
  integer:: iLower_D(nDim), iUpper_D(nDim) ! index limits for each box

  integer:: iLowerBc_D(nDim), iUpperBc_D(nDim) ! index limit for ghost cells
  integer:: jLowerBc_D(nDim), jUpperBc_D(nDim) ! index limit for phys. cells

  ! Mapping between boundaries: same index order and same direction
  integer, parameter:: iIndexMap_D(nDim) = (/0,1,2/) 
  integer, parameter:: iIndexDir_D(nDim) = (/1,1,1/)

  integer:: iStencil
  integer:: iStencil_I(nStencil) = (/ (iStencil, iStencil=0,nStencil-1) /)
  integer:: DiStencil_DI(nDim,nStencil)   ! Stencil description

  integer:: iObjectType

  integer          :: nValue     ! number of matrix elements for 1 box
  real, allocatable:: Value_I(:) ! matrix elements

  integer:: iVerboseAmg        = 0    ! 0..3
  integer:: MaxRowElementsAmg  = 0    ! 2, 4 or 6 (2*nDim)
  integer:: iCoarsenAmg        = 6    ! 0,1,3,6,7,8,9,10,11,21,22
  integer:: iRelaxAmg          = 6    ! 0..6,8,9,15..18
  integer:: iInterpolateAmg    = 0    ! 0..14
  real::    StrongThresholdAmg = 0.25 ! 0.25 for 2D, 0.5-0.6 for 3D
  real::    TruncFactorAmg     = 0.0  ! ?
  
contains
  !==========================================================================
  subroutine hypre_read_param

    use ModReadParam, ONLY: read_var

    call read_var('iVerboseAmg',        iVerboseAmg)
    call read_var('MaxRowElementsAmg',  MaxRowElementsAmg)
    call read_var('iCoarsenAmg',        iCoarsenAmg)
    call read_var('iRelaxAmg',          iRelaxAmg)
    call read_var('iInterpolateAmg',    iInterpolateAmg)  
    call read_var('StrongThresholdAmg', StrongThresholdAmg)
    call read_var('TruncFactorAmg',     TruncFactorAmg)   
    call read_var('UseSinglePart',      UseSinglePart)

  end subroutine hypre_read_param

  !==========================================================================
  subroutine hypre_initialize

    integer:: iValue, iR, iTheta, iPhi, iError

    character(len=*), parameter:: NameSub = 'hypre_initialize'
    !-------------------------------------------------------------------------

    if(DoTestMe)write(*,*) NameSub,' starting'

    if(NamePreconditioner == 'MG')then
       iObjectType = HYPRE_STRUCT
    else
       iObjectType = HYPRE_PARCSR
    end if

    if(UseSinglePart)then
       ! one part all together
       iPart = 0
       nPart = 1

       ! box = local domain with global index space
       iLower_D = (/  1, iTheta0+1,      iPhi0+1 /)
       iUpper_D = (/ nR, iTheta0+nTheta, iPhi0+nPhi /)
    else
       ! one part on each processor
       iPart = iProc
       nPart = nProc

       ! part = local domain with local index space
       iLower_D = (/  1, 1, 1 /)
       iUpper_D = (/ nR, nTheta, nPhi /)
    end if

    ! Create an empty 3D grid object
    call HYPRE_SStructGridCreate(iComm, nDim, nPart, i8Grid, iError)

    ! Add local box/part
    call HYPRE_SStructGridSetExtents(i8Grid, iPart, iLower_D, iUpper_D, iError)

    ! Single cell centered variable on all parts
    do jPart = 0, nPart - 1
       call HYPRE_SStructGridSetVariables(i8Grid, jPart, nVar, iVarType_I, iError)
    end do

    if(UseSinglePart)then

       ! This solution does not work, because not implemented for CSR storage.
       !nPeriod_D = (/0,0,nPhiAll/)
       !call HYPRE_SStructGridSetPeriodic(i8Grid, &
       !     iPart, nPeriod_D, iError )

       ! Setup periodic boundaries in Phi direction
       iLowerBc_D = (/  1,         1,    0 /)
       iUpperBc_D = (/ nR, nThetaAll,    0 /)
       jLowerBc_D = (/  1,         1, nPhiAll /)
       jUpperBc_D = (/ nR, nThetaAll, nPhiAll /)
       
       call HYPRE_SStructGridSetNeighborPart( i8Grid, &
            iPart, iLowerBc_D, iUpperBc_D, &
            iPart, jLowerBc_D, jUpperBc_D, &
            iIndexMap_D, iIndexDir_D, iError)
       
       iLowerBc_D = (/  1,         1, nPhiAll+1 /)
       iUpperBc_D = (/ nR, nThetaAll, nPhiAll+1 /)
       jLowerBc_D = (/  1,         1,         1 /)
       jUpperBc_D = (/ nR, nThetaAll,         1 /)
       
       call HYPRE_SStructGridSetNeighborPart( i8Grid, &
            iPart, iLowerBc_D, iUpperBc_D, &
            iPart, jLowerBc_D, jUpperBc_D, &
            iIndexMap_D, iIndexDir_D, iError)
    else

       ! Setup connection between parts

       ! Left neighbor
       if(iProcPhi > 0)then
          jPart = iProc - 1
       else
          jPart = iProc + nProcPhi - 1
       end if
       if(DoTestMe) write(*,*)'Left   iPart, jPart=',iPart, jPart

       iLowerBc_D = (/  1,      1,    0 /)
       iUpperBc_D = (/ nR, nTheta,    0 /)
       jLowerBc_D = (/  1,      1, nPhi /)
       jUpperBc_D = (/ nR, nTheta, nPhi /)

       call HYPRE_SStructGridSetNeighborPart( i8Grid, &
            iPart, iLowerBc_D, iUpperBc_D, &
            jPart, jLowerBc_D, jUpperBc_D, &
            iIndexMap_D, iIndexDir_D, iError)

       ! Right neighbor
       if(iProcPhi < nProcPhi - 1)then
          jPart = iProc + 1
       else
          jPart = iProc - (nProcPhi - 1)
       end if
       if(DoTestMe) write(*,*)'Right  iPart, jPart=', iPart, jPart

       iLowerBc_D = (/  1,      1, nPhi+1 /)
       iUpperBc_D = (/ nR, nTheta, nPhi+1 /)
       jLowerBc_D = (/  1,      1,      1 /)
       jUpperBc_D = (/ nR, nTheta,      1 /)

       call HYPRE_SStructGridSetNeighborPart( i8Grid, &
            iPart, iLowerBc_D, iUpperBc_D, &
            jPart, jLowerBc_D, jUpperBc_D, &
            iIndexMap_D, iIndexDir_D, iError)

       ! Bottom neighbor
       if(iProc >= nProcPhi)then
          jPart = iProc - nProcPhi
          if(DoTestMe) write(*,*)'Bottom iPart, jPart=', iPart, jPart

          iLowerBc_D = (/  1,      0,    1 /)
          iUpperBc_D = (/ nR,      0, nPhi /)
          jLowerBc_D = (/  1, nTheta,    1 /)
          jUpperBc_D = (/ nR, nTheta, nPhi /)

          call HYPRE_SStructGridSetNeighborPart( i8Grid, &
               iPart, iLowerBc_D, iUpperBc_D, &
               jPart, jLowerBc_D, jUpperBc_D, &
               iIndexMap_D, iIndexDir_D, iError)
       end if

       ! Top neighbor
       if(iProc < nProc - nProcPhi)then
          jPart = iProc + nProcPhi
          if(DoTestMe) write(*,*)'Top    iPart, jPart=', iPart, jPart

          iLowerBc_D = (/  1, nTheta+1,    1 /)
          iUpperBc_D = (/ nR, nTheta+1, nPhi /)
          jLowerBc_D = (/  1,        1,    1 /)
          jUpperBc_D = (/ nR,        1, nPhi /)

          call HYPRE_SStructGridSetNeighborPart( i8Grid, &
               iPart, iLowerBc_D, iUpperBc_D, &
               jPart, jLowerBc_D, jUpperBc_D, &
               iIndexMap_D, iIndexDir_D, iError)
       end if

    end if

    ! Assemble grid from all processors
    call HYPRE_SStructGridAssemble(i8Grid, iError)
    if(iError/=0)write(*,*)'ERROR: HYPRE_SStructGridAssemble failed'
    if(DoTestMe)write(*,*) NameSub,' HYPRE_SStructGridAssemble done'

    ! Define index offsets for the 7-point stencil
    DiStencil_DI(:,1) = (/  0, 0, 0 /)
    DiStencil_DI(:,2) = (/ -1, 0, 0 /)
    DiStencil_DI(:,3) = (/ +1, 0, 0 /)
    DiStencil_DI(:,4) = (/  0,-1, 0 /)
    DiStencil_DI(:,5) = (/  0,+1, 0 /)
    DiStencil_DI(:,6) = (/  0, 0,-1 /)
    DiStencil_DI(:,7) = (/  0, 0,+1 /)

    call HYPRE_SStructStencilCreate(nDim, nStencil, i8Stencil, iError)

    do iStencil = 1, nStencil
       call HYPRE_SStructStencilSetEntry(i8Stencil, &
            iStencil-1, DiStencil_DI(1,iStencil), iVar, iError)
    enddo

    ! Create the graph object
    call HYPRE_SStructGraphCreate(iComm, i8Grid, i8Graph, iError)

    ! Tell the graph which stencil to use for each variable on each part 
    do jPart = 0, nPart - 1
       call HYPRE_SStructGraphSetStencil(i8Graph, jPart, iVar, i8Stencil, iError)
    end do

    ! Assemble the graph
    call HYPRE_SStructGraphAssemble(i8Graph, iError)

    if(DoTestMe)write(*,*) NameSub,' HYPRE_SStructGraphAssemble done'

    ! Create an empty matrix object
    call HYPRE_SStructMatrixCreate(iComm, i8Graph, i8A, iError)

    ! Set storage type
    call HYPRE_SStructMatrixSetObjectTyp(i8A, iObjectType, iError)

    ! Get ready to set values
    call HYPRE_SStructMatrixInitialize(i8A, iError)
    if(DoTestMe)write(*,*) NameSub,' HYPRE_SStructMatrixInitialize done'

    ! Set non-zero matrix elements
    nValue = nStencil*nR*nTheta*nPhi
    allocate(Value_I(nValue))

    Value_I(1:nValue:nStencil) = -d_I
    Value_I(2:nValue:nStencil) = -e_I
    Value_I(3:nValue:nStencil) = -f_I
    Value_I(4:nValue:nStencil) = -e1_I
    Value_I(5:nValue:nStencil) = -f1_I
    Value_I(6:nValue:nStencil) = -e2_I
    Value_I(7:nValue:nStencil) = -f2_I

    call HYPRE_SStructMatrixSetBoxValues(i8A, iPart, iLower_D, iUpper_D, &
         iVar, nStencil, iStencil_I, Value_I, iError)

    ! Assemble matrix
    call HYPRE_SStructMatrixAssemble(i8A, iError)
    if(DoTestMe)write(*,*) NameSub,' HYPRE_SStructMatrixAssemble done'
    deallocate(Value_I)

    ! Create empty vector objects for RHS and solution
    call HYPRE_SStructVectorCreate(iComm, i8Grid, i8B, iError)
    call HYPRE_SStructVectorCreate(iComm, i8Grid, i8X, iError)

    ! Set object type for the vectors
    call HYPRE_SStructVectorSetObjectTyp(i8B, iObjectType, iError)
    call HYPRE_SStructVectorSetObjectTyp(i8X, iObjectType, iError)

    ! Initialize vectors
    call HYPRE_SStructVectorInitialize(i8B, iError)
    call HYPRE_SStructVectorInitialize(i8X, iError)

    ! Pass matrix to the solvers
    call HYPRE_SStructMatrixGetObject(i8A, i8ParA, iError)

    select case(NamePreconditioner)
    case('MG')
       ! Create the PFMG preconditioner
       call HYPRE_StructPFMGCreate(iComm, i8Precond, iError)

       ! Set PFMG parameters
       call HYPRE_StructPFMGSetMaxIter(i8Precond, 1, iError)
       call HYPRE_StructPFMGSetTol(i8Precond, 0.0, iError)
       call HYPRE_StructPFMGSetZeroGuess(i8Precond, iError)
       call HYPRE_StructPFMGSetNumPreRelax(i8Precond, 2, iError)
       call HYPRE_StructPFMGSetNumPostRelax(i8Precond, 2, iError)
       ! Non-Galerkin coarse grid (more efficient for this problem)
       call HYPRE_StructPFMGSetRAPType(i8Precond, 1, iError)
       ! R/B Gauss-Seidel
       call HYPRE_StructPFMGSetRelaxType(i8Precond, 2, iError)
       ! Skip relaxation on some levels (more efficient for this problem)
       call HYPRE_StructPFMGSetSkipRelax(i8Precond, 1, iError)

       if(DoTestMe)write(*,*) NameSub,' HYPRE_StructPFMGSetSkipRelax done'

    case('AMG')
       ! Create the BoomerAMG as a preconditioner
       call HYPRE_BoomerAMGCreate(i8Precond, iError)

       ! Set BoomerAMG parameters
       call HYPRE_BoomerAMGSetMaxIter(i8Precond, 1, iError)
       call HYPRE_BoomerAMGSetTol(i8Precond, 0.0, iError)
       call HYPRE_BoomerAMGSetNumSweeps(i8Precond, 1, iError)

       ! Adjustable parameters
       call HYPRE_BoomerAMGSetPrintLevel(   i8Precond, iVerboseAmg, iError)
       call HYPRE_BoomerAMGSetPMaxElmts(    i8Precond, MaxRowElementsAmg, iError)
       call HYPRE_BoomerAMGSetCoarsenType(  i8Precond, iCoarsenAmg, iError)
       call HYPRE_BoomerAMGSetRelaxType(    i8Precond, iRelaxAmg, iError)
       call HYPRE_BoomerAMGSetInterpType(   i8Precond, iInterpolateAmg, iError)
       call HYPRE_BoomerAMGSetStrongThrshld(i8Precond, StrongThresholdAmg, iError)
       call HYPRE_BoomerAMGSetTruncFactor(  i8Precond, TruncFactorAmg, iError)

       if(UsePreconditioner)then
          ! Setup AMG preconditioner for Krylov solver
          if(UseTiming)write(*,*)NameSub, &
               ' time before BoomerAMGSetup:', MPI_WTIME() - TimeStart
          call HYPRE_BoomerAMGSetup(i8Precond, i8ParA, i8ParB, i8ParX, iError)
          if(UseTiming)write(*,*)NameSub, &
               ' time after BoomerAMGSetup:', MPI_WTIME() - TimeStart
       end if

    end select

    if(DoTestMe)write(*,*) NameSub,' finished'

  end subroutine hypre_initialize

  !==========================================================================

  subroutine hypre_solver

    integer:: iValue, iR, iTheta, iPhi, iError

    character(len=*), parameter:: NameSub = 'hypre_solver'
    !-------------------------------------------------------------------------

    ! Set RHS values
    nValue = nR*nTheta*nPhi
    allocate(Value_I(nValue))

    if(DoTestMe)write(*,*) NameSub,' starting n, maxval, minval(Rhs_C)=', &
         nValue, maxval(Rhs_C), minval(Rhs_C), Rhs_C(1,1,1)

    iValue = 0
    do iPhi = 1, nPhi; do iTheta = 1, nTheta; do iR = 1, nR
       iValue = iValue + 1
       Value_I(iValue) = -Rhs_C(iR,iTheta,iPhi)
    end do; end do; end do

    call HYPRE_SStructVectorSetBoxValues(i8B, iPart, iLower_D, iUpper_D, &
         iVar, Value_I, iError)

    if(DoTestMe)write(*,*) NameSub,' set RHS'

    ! Set initial guess value to zero
    Value_I = 0.0
    call HYPRE_SStructVectorSetBoxValues(i8X, iPart, iLower_D, iUpper_D, &
         iVar, Value_I, iError)

    if(DoTestMe)write(*,*) NameSub,' set X=0'

    ! Assemble vectors
    call HYPRE_SStructVectorAssemble(i8X, iError)
    call HYPRE_SStructVectorAssemble(i8B, iError)
    if(DoTestMe)write(*,*) NameSub,' HYPRE_SStructVectorAssemble done'

    ! Pass vectors to the solvers
    call HYPRE_SStructVectorGetObject(i8B, i8ParB, iError)
    call HYPRE_SStructVectorGetObject(i8X, i8ParX, iError)
    if(DoTestMe)write(*,*) NameSub,' passed vectors to solver'

    select case(NameSolver)
    case('AMG')
       ! Create an empty BoomerAMG solver
       call HYPRE_BoomerAMGCreate(i8Solver, iError)
       if(DoTestMe)write(*,*) NameSub,' HYPRE_BoomerAMGCreate done'
       ! print solve info + parameters
       call HYPRE_BoomerAMGSetPrintLevel(i8Solver, 3, iError)
       ! Falgout coarsening
       call HYPRE_BoomerAMGSetCoarsenType(i8Solver, 6, iError)
       ! G-S/Jacobi hybrid relaxation
       call HYPRE_BoomerAMGSetRelaxType(i8Solver, 3, iError)
       ! Sweeeps on each level
       call HYPRE_BoomerAMGSetNumSweeps(i8Solver, 1, iError)
       ! maximum number of levels
       call HYPRE_BoomerAMGSetMaxLevels(i8Solver, 20, iError)
       ! conv. tolerance
       call HYPRE_BoomerAMGSetTol(i8solver, Tolerance, iError)
       if(DoTestMe)write(*,*) NameSub,' HYPRE_BoomerAMGSetTol done'

       ! Now setup and solve
       if(UseTiming)write(*,*)NameSub, &
            ' time before BoomerAMGSetup:', MPI_WTIME() - TimeStart

       call HYPRE_BoomerAMGSetup(i8Solver, i8ParA, i8ParB, i8ParX, iError)

       if(UseTiming)write(*,*)NameSub, &
            ' time before BoomerAMGSetup:', MPI_WTIME() - TimeStart

       call HYPRE_BoomerAMGSolve(i8Solver, i8ParA, i8ParB, i8ParX, iError)

       if(UseTiming)write(*,*)NameSub, &
            ' time after BoomerAMGSolve:', MPI_WTIME() - TimeStart

       ! Free memory
       call HYPRE_BoomerAMGDestroy(i8Solver, iError)

    case('GMRES')
       select case(NamePreconditioner)
       case('MG')
          ! Create an empty GMRES solver
          call HYPRE_StructGMRESCreate(iComm, i8Solver, iError)

          ! Set GMRES parameters
          call HYPRE_StructGMRESSetTol(i8Solver, Tolerance, iError)
          call HYPRE_StructGMRESSetPrintLevel(i8Solver, 2, iError) !!! 2
          call HYPRE_StructGMRESSetMaxIter(i8Solver, 50, iError)

          if(DoTestMe)write(*,*) NameSub,' HYPRE_StructGMRESSetMaxIter done'

          ! Set preconditioner (PFMG = 1) and solve
          call HYPRE_StructGMRESSetPrecond(i8Solver, 1, i8Precond, iError)
          if(DoTestMe)write(*,*) NameSub,' HYPRE_StructGMRESSetPrecond done'

          if(UseTiming)write(*,*)NameSub, &
               ' time before HYPRE_StructGMRESSetup:', MPI_WTIME() - TimeStart

          call HYPRE_StructGMRESSetup(i8Solver, i8ParA, i8ParB, i8ParX, iError)
          if(DoTestMe)write(*,*) NameSub,' HYPRE_StructGMRESSetup done'

          if(UseTiming)write(*,*)NameSub, &
               ' time after HYPRE_StructGMRESSetup:', MPI_WTIME() - TimeStart

          call HYPRE_StructGMRESSolve(i8Solver, i8ParA, i8ParB, i8ParX, iError)
          if(DoTestMe)write(*,*) NameSub,' HYPRE_StructGMRESSolve done'

          if(UseTiming)write(*,*)NameSub, &
               ' time after HYPRE_StructGMRESSolve:', MPI_WTIME() - TimeStart

          ! Free memory
          call HYPRE_StructGMRESDestroy(i8Solver, iError)
          call HYPRE_StructPFMGDestroy(i8Precond, iError)
       case('AMG')
          ! Create an empty GMRES solver
          call HYPRE_ParCSRGMRESCreate(iComm, i8Solver, iError)

          ! Set GMRES parameters
          call HYPRE_ParCSRGMRESSetTol(i8Solver, Tolerance, iError)
          call HYPRE_ParCSRGMRESSetPrintLevel(i8Solver, 100, iError) !!! 2
          call HYPRE_ParCSRGMRESSetMaxIter(i8Solver, 50, iError)

          if(DoTestMe)write(*,*) NameSub,' HYPRE_ParCSRGMRESSetMaxIter done'

          ! Set preconditioner (BoomerAMG = 2) and solve
          call HYPRE_ParCSRGMRESSetPrecond(i8Solver, 2, i8Precond, iError)
          if(DoTestMe)write(*,*) NameSub,' HYPRE_ParCSRGMRESSetPrecond done'

          if(UseTiming)write(*,*) NameSub, &
               ' time before HYPRE_ParCSRGMRESSetup:', MPI_WTIME() - TimeStart

          call HYPRE_ParCSRGMRESSetup(i8Solver, i8ParA, i8ParB, i8ParX, iError)
          if(DoTestMe)write(*,*) NameSub,' HYPRE_ParCSRGMRESSetup done'

          if(UseTiming)write(*,*)NameSub, &
               ' time after HYPRE_ParCSRGMRESSetup:', MPI_WTIME() - TimeStart

          call HYPRE_ParCSRGMRESSolve(i8Solver, i8ParA, i8ParB, i8ParX, iError)
          if(DoTestMe)write(*,*) NameSub,' HYPRE_ParCSRGMRESSolve done'

          if(UseTiming)write(*,*)NameSub, &
               ' time after HYPRE_ParCSRGMRESSolve:', MPI_WTIME() - TimeStart

          ! Free memory
          call HYPRE_ParCSRGMRESDestroy(i8Solver, iError)
          call HYPRE_BoomerAMGDestroy(i8Precond, iError)
       end select
    end select

    ! Get result
    Potential_C = 0.0
    call HYPRE_SStructVectorGather(i8x, iError);
    call HYPRE_SStructVectorGetBoxValues(i8X, iPart, &
         iLower_D, iUpper_D, iVar, Potential_C, iError)
    if(DoTestMe)write(*,*) NameSub,' HYPRE_SStructVectorGetBoxValues done'

    ! Free memory
    call HYPRE_SStructGridDestroy(i8Grid, iError)
    call HYPRE_SStructStencilDestroy(i8Stencil, iError)
    call HYPRE_SStructGraphDestroy(i8Graph, iError)
    call HYPRE_SStructMatrixDestroy(i8A, iError)
    call HYPRE_SStructVectorDestroy(i8B, iError)
    call HYPRE_SStructVectorDestroy(i8X, iError)

    deallocate(Value_I)

    if(DoTestMe)write(*,*) NameSub,' finished'

  end subroutine hypre_solver

  !===========================================================================

  subroutine hypre_preconditioner(n, y_I)

    integer, intent(in):: n
    real, intent(inout):: y_I(n)

    integer:: iError
    real, allocatable:: Value_I(:)

    logical, parameter:: DoDebug = .false.

    character(len=*), parameter:: NameSub = 'hypre_preconditioner'
    !-------------------------------------------------------------------------

    if(DoDebug)write(*,*) NameSub,' starting n, maxval, minval, y_I(1)=', &
         n, maxval(y_I), minval(y_I), y_I(1)

    ! Preconditioning: y'= AMG.y

    ! Set y_I as the RHS
    call HYPRE_SStructVectorSetBoxValues(i8B, iPart, iLower_D, iUpper_D, &
         iVar, y_I, iError)

    if(DoDebug)write(*,*) NameSub,' set RHS, iLower_D, iUpper_D=', &
         iLower_D, iUpper_D

    ! Set initial guess value to zero
    allocate(Value_I(n))
    Value_I = 0.0
    call HYPRE_SStructVectorSetBoxValues(i8X, iPart, iLower_D, iUpper_D, &
         iVar, Value_I, iError)
    deallocate(Value_I)

    if(DoDebug)write(*,*) NameSub,' set X=0'

    ! Assemble vectors
    call HYPRE_SStructVectorAssemble(i8X, iError)
    call HYPRE_SStructVectorAssemble(i8B, iError)

    if(DoDebug)write(*,*) NameSub,' HYPRE_SStructVectorAssemble done'

    ! Pass the vectors to the solvers
    call HYPRE_SStructVectorGetObject(i8B, i8ParB, iError)
    call HYPRE_SStructVectorGetObject(i8X, i8ParX, iError)

    if(DoDebug)write(*,*) NameSub,' passed vectors to AMG'

    call HYPRE_BoomerAMGSolve(i8Precond, i8ParA, i8ParB, i8ParX, iError)

    if(DoDebug)write(*,*) NameSub,' applied AMG preconditioner'

    ! Get back solution
    call HYPRE_SStructVectorGather(i8x, iError);
    call HYPRE_SStructVectorGetBoxValues(i8X, iPart, &
         iLower_D, iUpper_D, iVar, y_I, iError)

    if(DoDebug)write(*,*) NameSub,' finished'

  end subroutine hypre_preconditioner

end module ModHypre
!============================================================
subroutine read_hypre_param

  ! This is here to avoid circular dependencies

  use ModHypre, ONLY: hypre_read_param
  call hypre_read_param

end subroutine read_hypre_param
