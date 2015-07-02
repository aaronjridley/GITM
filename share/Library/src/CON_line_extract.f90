!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!BOP
!MODULE: CON_line_extract - extract field and stream lines in parallel
!INTERFACE:
module CON_line_extract

  !DESCRIPTION:
  ! The ray position and the MHD state along the rays (ray line) 
  ! can be collected into an array and sorted by the length coordinate. 
  ! This class provides the infrastructure for collecting, 
  ! sorting and providing the data extracted along multiple ray lines.

  !USES:
  use ModMpi
  use ModSort

  implicit none

  save    ! save all variables

  private ! except

  !PUBLIC MEMBER FUNCTIONS:
  public :: line_init         ! Initialize storage for ray lines
  public :: line_clean        ! Clean up storage for ray lines
  public :: line_put          ! Store state of a single point along a ray line
  public :: line_collect      ! Collect all ray line data onto 1 processor
  public :: line_get          ! Get all (sorted) ray line data from 1 processor
  public :: line_test         ! Unit tester

  !REVISION HISTORY:
  ! 09May04 - Gabor Toth <gtoth@umich.edu> - initial prototype/prolog/code
  !EOP

  ! Private constants
  character(len=*),  parameter :: NameMod='CON_line_extract'

  integer       :: nVar=0        ! Number of variables for each point
  integer       :: nPoint=0      ! Actual number of points
  integer       :: MaxPoint=0    ! Allocated number of points
  real, pointer :: State_VI(:,:) ! State at each point along all rays

contains

  !BOP ========================================================================
  !IROUTINE: line_init - initialize storage for ray lines
  !INTERFACE:
  subroutine line_init(nVarIn)

    !INPUT ARGUMENTS:
    integer, intent(in) :: nVarIn   ! Number of variables to store

    !DESCRIPTION:
    ! Initialize the ray line storage.
    !EOP

    character (len=*), parameter :: NameSub = NameMod//'::line_init'
    !-------------------------------------------------------------------------
    
    if(nVar >  0     ) call line_clean ! clean up previous allocation

    ! Set number of variables for each point
    nVar = nVarIn
    ! Initialize storaga
    nullify(State_VI)

  end subroutine line_init

  !BOP ========================================================================
  !IROUTINE: line_clean - clean storage for line data
  !INTERFACE:
  subroutine line_clean

    !DESCRIPTION:
    ! Clean the ray line storage.
    !EOP

    character (len=*), parameter :: NameSub = NameMod//'::line_clean'
    !-------------------------------------------------------------------------
    
    if(nVar == 0) RETURN  ! nothing to do, already clean

    nVar     = 0
    nPoint   = 0
    MaxPoint = 0
    if(associated(State_VI)) deallocate(State_VI)

  end subroutine line_clean

  !BOP ========================================================================
  !IROUTINE: line_extend - allocate or extend storage for line data
  !INTERFACE:
  subroutine line_extend(MaxPointIn)

    integer, intent(in) :: MaxPointIn

    real, pointer :: OldState_VI(:,:)
    !------------------------------------------------------------------------
    if(.not.associated(State_VI))then
       allocate(State_VI(0:nVar, MaxPointIn)) ! allocate storage
       MaxPoint = MaxPointIn                  ! set buffer size
    else
       OldState_VI => State_VI                ! store old values
       allocate(State_VI(0:nVar, MaxPointIn)) ! allocate new storage
       State_VI(:,1:nPoint) = &
            OldState_VI(:,1:nPoint)           ! copy old values
       deallocate(OldState_VI)                ! free old storage
       MaxPoint = MaxPointIn                  ! change buffer size
    end if
       
  end subroutine line_extend

  !BOP ========================================================================
  !IROUTINE: line_put - store state for a point along a ray line
  !INTERFACE:
  subroutine line_put(iLine, nVarIn, Line_V)

    !INPUT ARGUMENTS:
    integer, intent(in) :: iLine
    integer, intent(in) :: nVarIn
    real,    intent(in) :: Line_V(nVarIn)

    !DESCRIPTION:
    ! Store the ray index iLine and the state Line_V(1:nVarIn) 
    ! for the line iLine.
    !EOP

    character (len=*), parameter :: NameSub = NameMod//'line_put'
    !------------------------------------------------------------------------
    if(nVarIn /= nVar)then
       write(*,*) NameSub,' ERROR: nVarIn, nVar=',nVarIn, nVar
       call CON_stop(NameSub//' ERROR: incorrect number of variables!')
    end if

    if(nPoint >= MaxPoint) call line_extend(nPoint + 100)
    nPoint                     = nPoint + 1
    State_VI(0       , nPoint) = real(iLine)
    State_VI(1:nVarIn, nPoint) = Line_V

  end subroutine line_put

  !BOP ========================================================================
  !IROUTINE: line_collect - collect ray line data onto 1 processor
  !INTERFACE:
  subroutine line_collect(iComm, iProcTo)

    !INPUT ARGUMENTS:
    integer, intent(in) :: iComm      ! MPI communicator for the involved PE-s
    integer, intent(in) :: iProcTo    ! Rank of the receiving PE

    !DESCRIPTION:
    ! Collect line data from the sending processors to the receiving processor.
    ! The sending arrays are emptied, the receiving array is extended.
    ! The receiving PE may or may not contain data prior to the collect.
    !EOP

    integer, allocatable :: nPoint_P(:)

    integer, parameter :: iTag = 1
    character(len=*), parameter :: NameSub = NameMod//'line_collect'

    integer :: iProc, nProc, iProcFrom, Status_I(MPI_STATUS_SIZE), iError
    !-------------------------------------------------------------------------

    ! Set processor rank
    call MPI_COMM_RANK(iComm, iProc, iError)
    call MPI_COMM_SIZE(iComm, nProc, iError)

    ! Collect number of points onto the receiving processors
    
    allocate(nPoint_P(0:nProc-1))
    do iProcFrom = 0, nProc-1

       if(iProc == iProcFrom .and. iProc == iProcTo)then
          nPoint_P(iProcFrom) = nPoint
          CYCLE
       end if

       if(iProc == iProcFrom) call MPI_SEND(nPoint, &
            1, MPI_INTEGER, iProcTo, iTag, iComm, iError)

       if(iProc == iProcTo) call MPI_RECV(nPoint_P(iProcFrom), &
            1, MPI_INTEGER, iProcFrom, iTag, iComm, Status_I, iError)
    end do

    ! Extend buffer on receiving processor
    if(iProc == iProcTo) call line_extend(sum(nPoint_P) + 100)

    ! Receive ray line data
    do iProcFrom = 0, nProc-1
       if(iProc == iProcFrom .and. iProc == iProcTo) CYCLE

       ! Send ray line data
       if(iProc == iProcFrom .and. nPoint>0) then
          call MPI_SEND(State_VI(:,1:nPoint), &
               nPoint * (nVar+1), MPI_REAL, &
               iProcTo, iTag, iComm, iError)
          deallocate(State_VI)
          nPoint = 0
       end if

       ! Receive ray line data
       if(iProc == iProcTo .and. nPoint_P(iProcFrom)>0) then
          call MPI_RECV(State_VI(:,nPoint+1:nPoint + nPoint_P(iProcFrom)), &
               nPoint_P(iProcFrom) * (nVar+1), MPI_REAL, &
               iProcFrom, iTag, iComm, Status_I, iError)
          nPoint = nPoint + nPoint_P(iProcFrom)
       end if

    end do

    deallocate(nPoint_P)

  end subroutine line_collect

  !BOP ========================================================================
  !IROUTINE: line_get - obtain all the ray states
  !INTERFACE:
  subroutine line_get(nVarOut, nPointOut, Line_VI, DoSort)

    !INPUT/OUTPUT ARGUMENTS:
    integer, intent(inout)         :: nVarOut
    integer, intent(inout)         :: nPointOut

    !OUTPUT ARGUMENTS:
    real,    intent(out), optional :: Line_VI(0:nVarOut, nPointOut)

    !INPUT ARGUMENTS:
    logical, intent(in), optional  :: DoSort

    !DESCRIPTION:
    ! Obtain all the ray states stored. 
    !\begin{verbatim}
    ! First get the number of variables and points:
    !
    !     call line_get(nVarOut, nPointOut)
    !
    ! (Line_VI should not be present!), then 
    !
    !    allocate(Line_VI(0:nVarOut, nPointOut)
    ! 
    ! (the 0 index is needed for the line index) and then obtain the state
    !
    !    call line_get(nVarOut, nPointOut, Line_VI, DoSort = .true.)
    !
    !\end{verbatim}
    !EOP

    character (len=*), parameter :: NameSub = NameMod//'::line_get'

    integer, allocatable :: iSorted_I(:)
    real :: Factor
    !------------------------------------------------------------------------
    if(.not. present(Line_VI))then
       ! Provide size
       nVarOut   = nVar
       nPointOUt = nPoint
       RETURN
    end if

    ! Check size
    if(nVar /= nVarOut .or. nPoint /= nPointOut)then
       write(*,*) NameSub,' ERROR: nVar  , nVarOut  =', &
            nVar, nVarOut
       write(*,*) NameSub,' ERROR: nPoint, nPointOut=', &
            nPoint, nPointOut
       call CON_stop(NameSub//' ERROR: incorrect array size')
    end if
    if(present(DoSort))then
       if(DoSort)then
          allocate(iSorted_I(nPoint))
          ! sort by line index (0) and then by length along line (1)
          Factor = 0.9/maxval(State_VI(1,1:nPoint))
          call sort_quick(nPoint, &
               State_VI(0,1:nPoint) + Factor*State_VI(1,1:nPoint), iSorted_I)
          ! Apply sorting
          Line_VI(:,1:nPoint) = State_VI(:,iSorted_I)
          ! Copy back into State_VI so the sorting can be avoided next time.
          State_VI(:,1:nPoint)= Line_VI(:,1:nPoint)
          deallocate(iSorted_I)
       else
          Line_VI(:,1:nPoint) = State_VI(:,1:nPoint)
       end if
    else
       Line_VI(:,1:nPoint) = State_VI(:,1:nPoint)
    end if

    nPoint = 0
       
  end subroutine line_get

  !BOP ========================================================================
  !IROUTINE: line_test - unit tester
  !INTERFACE:
  subroutine line_test

    !DESCRIPTION:
    ! Test the CON_line_extract module. This subroutine should be called from
    ! a small stand alone program.
    !EOP

    integer, parameter :: nVarTest = 4
    integer :: iProc, nProc, iError
    integer :: nVarOut, nPointOut, iPoint, iProcFrom
    real    :: Length, Index, State_V(2:nVarTest)
    real, allocatable :: Result_VI(:,:)

    character(len=*), parameter :: NameSub=NameMod//'::line_test'
    !------------------------------------------------------------------------

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nProc, iError)
    call MPI_COMM_RANK(MPI_COMM_WORLD, iProc, iError)

    if(iProc==0)write(*,'(a)')'Testing line_init'
    call line_init(nVarTest)

    if(nVar/=4)write(*,*)NameSub,' line_init failed, iProc=',iProc, &
         ' nVar =',nVar,' should be 4 !'

    if(iProc==0)write(*,'(a)')'Testing line_put'

    call line_put(1, 4, (/2.0*iProc,   10.0+iProc, 20.0+iProc, 30.0+iProc/) )
    call line_put(1, 4, (/2.0*iProc+1, 40.0+iProc, 50.0+iProc, 60.0+iProc/) )
    call line_put(2, 4, (/2.0*iProc,   70.0+iProc, 80.0+iProc, 90.0+iProc/) )

    if(nPoint /= 3)write(*,*)NameSub,' line_put failed, iProc=',iProc, &
         ' nPoint =',nPoint,' should be 3'
    if(MaxPoint /= 100)write(*,*)NameSub,' line_put failed, iProc=',iProc, &
         ' MaxPoint =',MaxPoint,' should be 100'

    do iPoint = 1, nPoint
       if(any(State_VI(:,iPoint) /= (/ &
            (iPoint+1)/2 +0.0, &
            2.0*iProc + mod(iPoint+1,2), &
            30.0*iPoint-20.0+iProc, &
            30.0*iPoint-10.0+iProc, &
            30.0*iPoint+iProc/)) ) &
            write(*,*)NameSub,' line_put failed, iProc=',iProc, &
            ' iPoint=',iPoint,' State_VI=',State_VI(:,iPoint)
    end do

    if(iProc==0)write(*,'(a)')'Testing line_collect'
    call line_collect(MPI_COMM_WORLD,0)

    if(iProc == 0)then
       if(nPoint /= nProc*3) &
            write(*,*)NameSub,' line_collect failed, iProc=',iProc, &
            ' nPoint =',nPoint,' should be',nProc*3
       if(MaxPoint /= 100+nProc*3) &
            write(*,*)NameSub,' line_collect failed, iProc=',iProc, &
            ' MaxPoint =',MaxPoint,' should be',100+nProc*3
    else
       if(nPoint /= 0) &
            write(*,*)NameSub,' line_collect failed, iProc=',iProc, &
            ' nPoint =',nPoint,' should be 0'
       if(MaxPoint /= 100) &
            write(*,*)NameSub,' line_collect failed, iProc=',iProc, &
            ' MaxPoint =',MaxPoint,' should be 100'

    end if

    if(iProc==0)then
       write(*,'(a)')'Testing line_get on PE 0'

       call line_get(nVarOut, nPointOut)
       if(nVarOut /= nVar)write(*,*)NameSub,' line_get failed', &
            ' nVarOut =',nVarOut,' should be ',nVar

       if(nPointOut /= nPoint)write(*,*)NameSub,' line_get failed', &
            ' nPointOut =',nPointOut,' should be', nPoint

       allocate(Result_VI(0:nVarOut,nPointOUt))
       call line_get(nVarOut, nPointOut, Result_VI, .true.)

       do iPoint = 1, nPointOut
          ! Test index
          if(iPoint <= 2*nProc)then
             Index = 1.0
          else
             Index = 2.0
          end if
          if(Result_VI(0,iPoint) /= Index) &
               write(*,*)NameSub,' line_get failed at iPoint=',iPoint, &
               ' Index=',Result_VI(0,iPoint),' should be ',Index

          ! Test length
          if(Index == 1.0)then
             Length = iPoint - 1.0
          else
             Length = 2.0*(iPoint - 2*nProc - 1)
          end if
          if(Result_VI(1,iPoint) /= real(Length)) &
               write(*,*)NameSub,' line_get failed at iPoint=',iPoint, &
               ' Length=',Result_VI(1,iPoint),' should be ',Length

          ! Test rest of the variables
          if(Index == 1.0)then
             iProcFrom = (iPoint - 1)/2
             if(mod(iPoint,2)==1)then
                State_V = (/10.0+iProcFrom, 20.0+iProcFrom, 30.0+iProcFrom/)
             else
                State_V = (/40.0+iProcFrom, 50.0+iProcFrom, 60.0+iProcFrom/)
             end if
          else
             iProcFrom = nint(Length / 2.0)
             State_V =    (/70.0+iProcFrom, 80.0+iProcFrom, 90.0+iProcFrom/)
          end if
          if(any(Result_VI(2:nVar,iPoint) /= State_V)) &
               write(*,*)NameSub,' line_get failed at iPoint=',iPoint, &
               ' State = ',Result_VI(2:nVar,iPoint),' should be ',State_V

       end do
    end if

    if(iProc==0)write(*,'(a)')'Testing line_clean'
    call line_clean

    if( nPoint/=0 .or. MaxPoint/=0 .or. associated(State_VI)) &
         write(*,*) NameSub,' line_clean failed, iProc=',iProc, &
         ' nPoint =',nPoint,' and Maxpoint =',MaxPoint,' should be 0', &
         ' associated(State_VI)=',associated(State_VI),' should be False!'
    
  end subroutine line_test

  !===========================================================================

end module CON_line_extract
