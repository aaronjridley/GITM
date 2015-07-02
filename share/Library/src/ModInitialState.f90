!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModInitialState

  ! Initialize a 2D domain with a finite number of material states.
  ! An arbitrary subset of the state variables can be set for each state.
  ! The rest of the variables are set to zero.
  !
  ! The material states are separated with an arbitrary number of 
  ! straight segments. For each segment the name of the material
  ! state is given for both sides.
  ! 
  ! The code can calculate the state for an arbitrary location.
  ! It can also set the levelset functions based on the distance
  ! to the closest material interface. Note that interfaces
  ! between different states but identical material have no effect
  ! on the levelset function. 

  implicit none

  SAVE

  private ! except

  public:: init_initial_state        ! allocate arrays
  public:: clean_initial_state       ! deallocate arrays
  public:: read_initial_state_param  ! read parameters
  public:: get_initial_state         ! set initial state at a point 
  public:: test_initial_state        ! unit test

  ! Local variables

  integer:: nVar = 0          ! Number of variables representing a state
  character(len=20), allocatable:: NameVar_V(:) ! State variable names
  integer:: iVarMaterial0         ! Variable index of first material
  integer:: nMaterial             ! Number of different materials

  integer:: nSegment = 0
  real,    allocatable:: LengthSegment_I(:)
  real,    allocatable:: StartSegment_DI(:,:)
  real,    allocatable:: VectorSegment_DI(:,:)
  integer, allocatable:: iLevelSegment_SI(:,:)
  integer, allocatable:: iStateSegment_SI(:,:)

  integer:: nMaterialState = 0 ! number of different material states
  character(len=10), allocatable:: NameMaterialState_I(:)
  real, allocatable:: MaterialState_VI(:,:)  ! state values

contains

  !===========================================================================

  subroutine init_initial_state(NameVarState, iVarMaterial1In, nMaterialIn)

    use ModUtilities, ONLY: split_string, lower_case

    character(len=*),  intent(in):: NameVarState

    integer, optional, intent(in):: iVarMaterial1In
    integer, optional, intent(in):: nMaterialIn

    integer:: iVar
    integer, parameter:: MaxVar = 100
    character(len=20), allocatable:: String_I(:)

    character(len=*), parameter:: NameSub='init_initial_state'
    !-------------------------------------------------------------------------
    ! Check if the module has been initialized already
    if(nVar > 0) RETURN

    allocate(String_I(MaxVar))
    call split_string(NameVarState, MaxVar, String_I, nVar)

    if(nVar < 1) call CON_stop(NameSub//': nVar < 1')

    allocate(NameVar_V(nVar))
    NameVar_V = String_I(1:nVar)
    deallocate(String_I)

    do iVar = 1, nVar
       call lower_case(NameVar_V(iVar))
    end do

    iVarMaterial0 = -1
    if(present(iVarMaterial1In)) iVarMaterial0 = iVarMaterial1In - 1

    nMaterial = 0
    if(present(nMaterialIn)) nMaterial = nMaterialIn

    if(nMaterial > 0 .and. &
         (iVarMaterial0 < 0 .or. iVarMaterial0 + nMaterial > nVar)) &
         call CON_stop(NameSub//': iVarMaterial0 has impossible value')

  end subroutine init_initial_state

  !===========================================================================

  subroutine clean_initial_state

    if(allocated(NameVar_V)) deallocate( &
         NameVar_V)

    nVar          = 0
    iVarMaterial0 = -1
    nMaterial     =  0

    if(allocated(LengthSegment_I)) deallocate( &
         LengthSegment_I,  &
         StartSegment_DI,  &
         VectorSegment_DI, &
         iLevelSegment_SI, &
         iStateSegment_SI)

    nSegment = 0

    if(allocated(NameMaterialState_I)) deallocate( &
         NameMaterialState_I, &
         MaterialState_VI)

    nMaterialState = 0

  end subroutine clean_initial_state

  !===========================================================================

  subroutine read_initial_state_param(NameCommand)

    use ModReadParam, ONLY: read_var
    use ModUtilities, ONLY: split_string

    character(len=*), intent(in):: NameCommand

    character(len=100):: String
    real:: Start_D(2), End_D(2)
    character(len=20):: NameMaterial1, NameMaterial2

    character(len=20), allocatable:: NameStateVar_I(:)
    integer, allocatable:: iVar_I(:)

    integer:: i, iVar, nStateVar, iError

    character(len=*), parameter:: NameSub='read_initial_state_param'
    !------------------------------------------------------------------------
    select case(NameCommand)
    case("#STATEDEFINITION")
       allocate(NameStateVar_I(nVar))

       ! Names of non-zero state variables
       call read_var("StringStateVar", String, &
            IsLowerCase=.true., DoReadWholeLine=.true.)

       ! Remove parameter name "string...." if present.
       i = index(String, 'string')
       if(i > 1)String(i:100) = ' '
       call split_string(String, nVar, NameStateVar_I, nStateVar)

       ! Find indexes for corresponding state variables
       allocate(iVar_I(nStateVar))
       LOOPSTATE: do i = 1, nStateVar
          do iVar = 1, nVar
             if(NameStateVar_I(i) == NameVar_V(iVar))then
                iVar_I(i) = iVar
                CYCLE LOOPSTATE
             end if
          end do
          call CON_stop(NameSub//': could not match variable name='// &
               trim(NameStateVar_I(i)))
       end do LOOPSTATE

       ! Number of different states in the initial condition
       call read_var('nMaterialState', nMaterialState) 

       ! Read the values for each material state
       allocate(NameMaterialState_I(nMaterialState))
       allocate(MaterialState_VI(nVar,nMaterialState))
       do i = 1, nMaterialState
          MaterialState_VI(:,i) = 0.0
          call read_var("Name State", String, &
               IsLowerCase=.true., DoReadWholeLine=.true.)
          read(String,*,IOSTAT=iError) &
               NameMaterialState_I(i), MaterialState_VI(iVar_I,i)
          if(iError /= 0)call CON_stop(NameSub// &
               ' could not read Name, State from '//trim(String))
       end do
       deallocate(NameStateVar_I, iVar_I)

    case("#STATEINTERFACE")
       call read_var('nSegment', nSegment)     ! Number of segments
       if(nSegment < 1) RETURN
       allocate(                            &
            LengthSegment_I(nSegment),    &
            StartSegment_DI(2,nSegment),  &
            VectorSegment_DI(2,nSegment), &
            iLevelSegment_SI(2,nSegment), &
            iStateSegment_SI(2,nSegment))

       do i = 1, nSegment
          call read_var("State1 State2 x1 y1 x2 y2", String, &
               IsLowerCase=.true., DoReadWholeLine=.true.)
          read(String,*,IOSTAT=iError) &
               NameMaterial1, NameMaterial2, Start_D, End_D
          if(iError /= 0)call CON_stop(NameSub// &
               ' could not read Material1 Material2 x1 y1 x2 y2 from '// &
               trim(String))

          iStateSegment_SI(1,i) = i_state(NameMaterial1)
          iStateSegment_SI(2,i) = i_state(NameMaterial2)

          if(nMaterial > 0)then
             iLevelSegment_SI(1,i) = i_level(NameMaterial1(1:2))
             iLevelSegment_SI(2,i) = i_level(NameMaterial2(1:2))
          end if

          StartSegment_DI(:,i)  = Start_D
          VectorSegment_DI(:,i) = End_D - Start_D
          LengthSegment_I(i)    = sqrt(sum( (End_D - Start_D)**2 ))
       end do
    case default
       call CON_stop(NameSub//': unknown command='//NameCommand)
    end select

  contains
    !=========================================================================
    integer function i_level(NameMaterial)

      character(len=2), intent(in):: NameMaterial
      integer:: iMaterial, iVar
      !-----------------------------------------------------------------------
      do iMaterial = 1, nMaterial
         iVar = iVarMaterial0 + iMaterial
         if(NameMaterial == NameVar_V(iVar)) then
            i_level = iMaterial
            RETURN
         end if
      end do

      call CON_stop(NameSub// &
           ' could not match material name='//NameMaterial)

    end function i_level

    !=========================================================================
    integer function i_state(NameMaterial)

      character(len=*), intent(in):: NameMaterial
      integer:: i
      !-----------------------------------------------------------------------
      if(nMaterialState == 0) call CON_stop(NameSub// &
           ' material states are not defined')

      do i = 1, nMaterialState
         if(NameMaterial == NameMaterialState_I(i)) then
            i_state = i
            RETURN
         end if
      end do

      write(*,*) NameSub,': nMaterialState      = ', nMaterialState
      write(*,*) NameSub,': NameMaterialState_I = ', NameMaterialState_I

      call CON_stop(NameSub// &
           ' could not match material state name='//NameMaterial)

    end function i_state

  end subroutine read_initial_state_param

  !===========================================================================

  subroutine get_initial_state(CoordIn_D, State_V, DoTest)

    real,    intent(in) :: CoordIn_D(2)          ! Coordinates of point
    real,    intent(out):: State_V(nVar)   ! State at this position

    logical, optional, intent(in):: DoTest

    ! Calculate levelset functions as the distance to the closest
    ! interface for all materials for a point located at CoordIn_D.

    integer:: iSegment      ! segment index
    integer:: iSide         ! which side of the segment
    integer:: iState        ! index of state the point belongs to
    integer:: iLevel        ! Level index on the side of the point
    integer:: jLevel        ! Level index on the opposite side
    real:: Coord_D(2), Segment_D(2), Point_D(2)
    real:: Length, Projection
    real:: d                ! Distance between segment and point
    real:: dNormal          ! Normal distance square between line and point
    real:: dMin             ! Distance to closest segment
    real, allocatable:: dMin_I(:) ! Distance to materials

    character(len=*), parameter:: NameSub = 'get_initial_state'
    !-------------------------------------------------------------------------

    ! Initialize state
    iState = 1
    State_V = 0.0
    dMin = huge(1.0)
    if(nMaterial > 0) then
       allocate(dMin_I(nMaterial))
       dMin_I = -huge(1.0)
    end if

    ! Check distance from all the material interface segments
    do iSegment = 1, nSegment

       ! point coordinate relative to the starting point
       Coord_D = CoordIn_D - StartSegment_DI(:,iSegment)

       ! Segment vector
       Segment_D = VectorSegment_DI(:,iSegment)

       ! Length of the segment
       Length = LengthSegment_I(iSegment)

       ! Distance of the point projected onto the line from point 1
       Projection = sum( Segment_D*Coord_D ) / Length

       ! Normal distance
       dNormal = sqrt(max(0.0,sum(Coord_D**2) - Projection**2))

       ! Find projected point relative to Point1 limited by the segments
       Point_D = min(1.0, max(0.0, Projection/Length))*Segment_D

       ! Distance to the projected point
       d = sqrt( sum((Coord_D - Point_D)**2) )

       ! Avoid problems around segments meeting at acute angles:
       ! when two segments are at the same distance from the point
       ! because their end points coincide, the one with the LARGER
       ! normal distance should be used. We slightly modify the distance
       ! (only when d is different from dNormal) to achieve this
       d = d + 1e-6*(d-dNormal)

       ! use cross product to figure out which side of the segment the point is
       iSide = 1
       if(Coord_D(1)*Segment_D(2) - Coord_D(2)*Segment_D(1) >= 0.0) iSide = 2

       if(d < dMin)then
          dMin = d
          iState = iStateSegment_SI(iSide,iSegment)
       end if

       if(nMaterial > 0)then
          ! Get the levels on the two sides of the segment.

          iLevel = iLevelSegment_SI(iSide,  iSegment)
          jLevel = iLevelSegment_SI(3-iSide,iSegment)

          if(iLevel /= jLevel)then
             ! Set level set function if distance is smaller than earlier
             if(d < abs(dMin_I(iLevel))) dMin_I(iLevel) =  d
             if(d < abs(dMin_I(jLevel))) dMin_I(jLevel) = -d
          end if
       end if

       if(.not.present(DoTest)) CYCLE
       if(.not.DoTest)          CYCLE

       write(*,*)NameSub,' iSegment, CoordIn_D=', iSegment, CoordIn_D
       write(*,*)NameSub,' Length,StartSegment=', &
            Length, StartSegment_DI(:,iSegment)
       write(*,*)NameSub,' Segment_D,Coord_D  =', Segment_D, Coord_D
       write(*,*)NameSub,' Projection, Point_D=', Projection, Point_D
       write(*,*)NameSub,' d, iLevel, jLevel  =', d, iLevel, jLevel

    end do

    State_V = MaterialState_VI(:,iState)
    if(nMaterial > 0)then
       State_V(iVarMaterial0+1:iVarMaterial0+nMaterial) = dMin_I
       deallocate(dMin_I)
    end if

  end subroutine get_initial_state

  !=========================================================================

  subroutine test_initial_state

    use ModPlotFile, ONLY: save_plot_file
    use ModReadParam, ONLY: read_file, read_init, read_line, read_command, &
         read_echo_set
    use ModMpi

    integer:: iProc, iError

    integer, parameter:: nI = 100, nJ = 50
    real, allocatable:: State_VC(:,:,:)
    real :: Coord_D(2)
    integer:: i, j

    character(len=100):: NameCommand

    character(len=*), parameter:: NameVar = &
         'Rho Xe Be Pl Au Ay Ux Uy Uz p level'
    integer, parameter:: iMaterial1 = 2, nMaterialTest = 5
    integer, parameter:: iMaterialLast = iMaterial1-1+nMaterialTest

    character(len=*), parameter:: NameSub = 'test_initial_state'
    !----------------------------------------------------------------------
    call MPI_COMM_RANK(MPI_COMM_WORLD, iProc, iError)

    call init_initial_state(NameVar, iMaterial1, nMaterialTest)

    ! Now nVar is set based on the list of variables
    allocate(State_VC(nVar,nI,nJ))

    ! Read in the parameters describing the material interfaces
    call read_file('test_initial_state.in')
    call read_init
    if(iProc==0)call read_echo_set(.true.)
    do
       if(.not.read_line()) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case("#STATEDEFINITION", "#STATEINTERFACE")
          call read_initial_state_param(NameCommand)
       case default
          call CON_stop(NameSub//': unknown command='//NameCommand)
       end select
    end do

    ! Set the levelset values for all cell centers
    do j = 1, nJ; do i = 1, nI
       ! Take a simple uniform grid (does not matter really)
       Coord_D = (/ i-0.5, j-0.5 /)

       call get_initial_state(Coord_D, State_VC(:,i,j))

       if(nMaterialTest < 1) CYCLE

       ! Limit smallest value (for materials that do not occur at all)
       State_VC(iMaterial1:iMaterialLast,i,j) = max(-777.7, &
            State_VC(iMaterial1:iMaterialLast,i,j))

       ! Set the level function based on the maximum level
       State_VC(nVar,i,j) = maxloc( &
            State_VC(iMaterial1:iMaterialLast,i,j), DIM=1)

    end do; end do

    if(iProc==0) call save_plot_file('test_initial_state.out', &
         NameVarIn = "x y "//NameVar//" level nSegment", &
         ParamIn_I = (/ nSegment+0.0 /), &
         CoordMinIn_D = (/ 0.5, 0.5 /), &
         CoordMaxIn_D = (/ nI-0.5, nJ-0.5 /), &
         VarIn_VII    = State_VC)

    deallocate(State_VC)

  end subroutine test_initial_state

end module ModInitialState
