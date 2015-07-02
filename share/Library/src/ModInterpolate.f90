!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModInterpolate

  ! Calculate second order accurate interpolation for 
  !
  ! - a uniform grid with normalized coordinates, or 
  ! - non-uniform grid with actual coordinates, or
  ! - any mixture of the two, i.e. only some of the coordinates are uniform
  !
  ! Normalized coordinates mean that the coordinates coincide with the 
  ! indexes at the grid points. For uniform grid this is a very fast algorithm.
  ! For non-uniform grid a binary search is needed. The coordinates are assumed
  ! to be either monotone increasing or monotone decreasing. 
  !
  ! One can interpolate both scalar and vector valued arrays.
  !
  ! If the coordinates are outside the allowed ranges and the DoExtrapolate
  ! argument is not present the code stops. If the DoExtrapolate argument
  ! is present and false, the last grid cell value is used. If DoExtrapolate
  ! is present and true, second order extrapolation is used.

  ! Examples of usage:
  !
  ! Cell based 2D uniform grid with ghost cells, scalar valued:
  !
  !     InterpolatedValue = bilinear(Value_II, 0, nI+1, 0, nJ+1, &
  !                         (/ (x - x0)/DeltaX, (y - y0)/DeltaY) /) )
  !
  ! Node based 2D grid with x(1)=y(1)=0.0, vector valued:
  !
  !     InterpolatedValue_V = bilinear(Value_VII, nVar, 1, nI, 1, nJ, &
  !                        (/ x/DeltaX+1, y/DeltaY+1 /) )
  !
  ! Nonuniform 3D grid with ghost cells, third coordinate is uniform, 
  ! scalar valued:
  !
  !     InterpolatedValue = trilinear(Value_III, -1, nI+2, -1, nJ+2, -1, nK+2,&
  !                       (/ x, y, (z - z0)/DeltaZ /), x_I, y_I)
  !

  implicit none

  private ! except

  public :: linear             ! 2nd order interpolation in 1D
  public :: bilinear           ! 2nd order interpolation in 2D
  public :: trilinear          ! 2nd order interpolation in 3D
  public :: find_cell          ! find cell in non-uniform grid
  public :: fit_parabola       ! fit a parabola around an extremum
  public :: test_interpolation ! unit test

  character(len=*), parameter :: NameMod='ModInterpolate'

  interface linear
     module procedure linear_scalar, linear_vector
  end interface

  interface bilinear
     module procedure bilinear_scalar, bilinear_vector
  end interface

  interface trilinear
     module procedure trilinear_scalar, trilinear_vector
  end interface

contains

  !=========================================================================
  real function linear_scalar(a_I, iMin, iMax, x, x_I, DoExtrapolate, &
       iCell, Dist)

    ! Calculate linear interpolation of a_I at position x
    ! Assume normalized coordinates unless x_I is present.
    ! If present x_I contains the coordinates in an increasing order.

    integer, intent(in) :: iMin, iMax
    real, intent(in)    :: a_I(iMin:iMax)
    real, intent(in)    :: x

    real,    intent(in), optional :: x_I(iMin:iMax)
    logical, intent(in), optional :: DoExtrapolate
    integer, intent(in), optional :: iCell
    real,    intent(in), optional :: Dist

    integer :: i1, i2
    real    :: Dx1, Dx2
    character (len=*), parameter :: NameSub=NameMod//'::linear_scalar'
    !--------------------------------------------------------------------------

    if ( present(iCell) .and. present(Dist) ) then

       ! Calculate the remaining cell indices and interpolation weights

       i1 = iCell
       i2 = i1 + 1
       Dx1 = Dist
       Dx2 = 1.0 - Dx1

    else

       ! Locate the point Xyz_D on the grid
       call find_cell(iMin, iMax, x, i1, Dx1, x_I, DoExtrapolate, &
            "Called from "//NameSub)

       ! Calculate the remaining cell indices and interpolation weights
       i2 = i1 + 1; Dx2 = 1.0 - Dx1

    end if

    ! Perform interpolation (or extrapolation)
    linear_scalar = (Dx2)*a_I(i1) + Dx1*a_I(i2)

  end function linear_scalar

  !=========================================================================
  function linear_vector(a_VI, nVar, iMin, iMax, x, x_I, DoExtrapolate, &
       iCell, Dist)

    ! Calculate linear interpolation of a_VI at position x
    ! Assume normalized coordinates unless x_I is present.
    ! If present x_I contains the coordinates in an increasing order.

    integer, intent(in) :: nVar, iMin, iMax
    real, intent(in)    :: a_VI(nVar, iMin:iMax)
    real, intent(in)    :: x

    real,    intent(in), optional :: x_I(iMin:iMax)
    logical, intent(in), optional :: DoExtrapolate
    integer, intent(in), optional :: iCell
    real,    intent(in), optional :: Dist

    ! return value
    real                :: linear_vector(nVar)

    integer :: i1, i2
    real    :: Dx1, Dx2
    character (len=*), parameter :: NameSub=NameMod//'::linear_vector'
    !--------------------------------------------------------------------------

    if ( present(iCell) .and. present(Dist) ) then

       ! Calculate the remaining cell indices and interpolation weights

       i1 = iCell
       i2 = i1 + 1
       Dx1 = Dist
       Dx2 = 1.0 - Dx1

    else

       ! Locate the point Xyz_D on the grid
       call find_cell(iMin, iMax, x, i1, Dx1, x_I, DoExtrapolate, &
            "Called from "//NameSub)

       ! Calculate the remaining cell indices and interpolation weights
       i2 = i1 + 1; Dx2 = 1.0 - Dx1

    end if

    ! Perform interpolation (or extrapolation) for multiple variables
    linear_vector = (Dx2)*a_VI(:,i1) + Dx1*a_VI(:,i2)

  end function linear_vector

  !=========================================================================
  real function bilinear_scalar( &
       a_II, iMin, iMax, jMin, jMax, Xy_D, x_I, y_I, DoExtrapolate, iCell_D, &
       Dist_D)

    ! Calculate bilinear interpolation of a_II at position Xy_D
    ! Assume normalized coordinates unless x_I and/or y_I are present.
    ! If present x_I and y_I contain the coordinates in an increasing order.

    integer, intent(in) :: iMin, iMax, jMin, jMax
    real, intent(in)    :: a_II(iMin:iMax,jMin:jMax)
    real, intent(in)    :: Xy_D(2)

    real,    intent(in), optional :: x_I(iMin:iMax)
    real,    intent(in), optional :: y_I(jMin:jMax)
    logical, intent(in), optional :: DoExtrapolate
    integer, intent(in), optional :: iCell_D(2)
    real,    intent(in), optional :: Dist_D(2)

    integer :: i1, i2, j1, j2
    real    :: Dx1, Dx2, Dy1, Dy2
    character (len=*), parameter :: NameSub=NameMod//'::bilinear_scalar'
    !--------------------------------------------------------------------------

    if ( present(iCell_D) .and. present(Dist_D) ) then

       ! Calculate the remaining cell indices and interpolation weights

       i1 = iCell_D(1)
       i2 = i1 + 1
       Dx1 = Dist_D(1)
       Dx2 = 1.0 - Dx1

       j1 = iCell_D(2)
       j2 = j1 + 1
       Dy1 = Dist_D(2)
       Dy2 = 1.0 - Dy1

    else

       ! Locate the point Xyz_D on the grid
       call find_cell(iMin, iMax, Xy_D(1), i1, Dx1, x_I, DoExtrapolate, &
            "Called for coord1 from "//NameSub)

       call find_cell(jMin, jMax, Xy_D(2), j1, Dy1, y_I, DoExtrapolate, &
            "Called for coord2 from "//NameSub)

       ! Calculate the remaining cell indices and interpolation weights
       i2 = i1 + 1; Dx2 = 1.0 - Dx1
       j2 = j1 + 1; Dy2 = 1.0 - Dy1

    end if

    ! Perform interpolation (or extrapolation)
    bilinear_scalar = Dy2*( Dx2*a_II(i1,j1)   &
         +                  Dx1*a_II(i2,j1))  &
         +            Dy1*( Dx2*a_II(i1,j2)   &
         +                  Dx1*a_II(i2,j2))

  end function bilinear_scalar

  !=========================================================================
  function bilinear_vector( &
       a_VII, nVar, iMin, iMax, jMin, jMax, Xy_D, x_I, y_I, DoExtrapolate, &
       iCell_D, Dist_D)

    ! Calculate bilinear interpolation of a_VII at position Xy_D
    ! Assume normalized coordinates unless x_I and/or y_I are present.
    ! If present x_I and y_I contain the coordinates in an increasing order.

    integer, intent(in) :: nVar, iMin, iMax, jMin, jMax
    real, intent(in)    :: a_VII(nVar, iMin:iMax,jMin:jMax)
    real, intent(in)    :: Xy_D(2)

    real,    intent(in), optional :: x_I(iMin:iMax)
    real,    intent(in), optional :: y_I(jMin:jMax)
    logical, intent(in), optional :: DoExtrapolate
    integer, intent(in), optional :: iCell_D(2)
    real,    intent(in), optional :: Dist_D(2)

    ! return value
    real                :: bilinear_vector(nVar)

    integer :: i1, i2, j1, j2
    real :: Dx1, Dx2, Dy1, Dy2
    character (len=*), parameter :: NameSub=NameMod//'::bilinear_vector'
    !--------------------------------------------------------------------------

    if ( present(iCell_D) .and. present(Dist_D) ) then

       ! Calculate the remaining cell indices and interpolation weights

       i1 = iCell_D(1)
       i2 = i1 + 1
       Dx1 = Dist_D(1)
       Dx2 = 1.0 - Dx1

       j1 = iCell_D(2)
       j2 = j1 + 1
       Dy1 = Dist_D(2)
       Dy2 = 1.0 - Dy1

    else

       ! Locate the point Xyz_D on the grid
       call find_cell(iMin, iMax, Xy_D(1), i1, Dx1, x_I, DoExtrapolate, &
            "Called for coord1 from "//NameSub)

       call find_cell(jMin, jMax, Xy_D(2), j1, Dy1, y_I, DoExtrapolate, &
            "Called for coord2 from "//NameSub)

       ! Calculate the remaining cell indices and interpolation weights
       i2 = i1 + 1; Dx2 = 1.0 - Dx1
       j2 = j1 + 1; Dy2 = 1.0 - Dy1

    end if

    ! Perform interpolation (or extrapolation) for multiple variables
    bilinear_vector = Dy2*( Dx2*a_VII(:,i1,j1)   &
         +                  Dx1*a_VII(:,i2,j1))  &
         +            Dy1*( Dx2*a_VII(:,i1,j2)   &
         +                  Dx1*a_VII(:,i2,j2))

  end function bilinear_vector

  !=========================================================================
  real function trilinear_scalar( &
       a_III, iMin, iMax, jMin, jMax, kMin, kMax, Xyz_D, &
       x_I, y_I, z_I, DoExtrapolate, iCell_D, Dist_D)

    ! Calculate trilinear interpolation of a_III at position Xyz_D

    integer, intent(in) :: iMin, iMax, jMin, jMax, kMin, kMax
    real, intent(in)    :: a_III(iMin:iMax,jMin:jMax,kMin:kMax)
    real, intent(in)    :: Xyz_D(3)

    real,    intent(in), optional :: x_I(iMin:iMax)
    real,    intent(in), optional :: y_I(jMin:jMax)
    real,    intent(in), optional :: z_I(kMin:kMax)
    logical, intent(in), optional :: DoExtrapolate

    integer,    intent(in), optional :: iCell_D(3)
    real,    intent(in), optional :: Dist_D(3)

    integer :: i1, i2, j1, j2, k1, k2
    real    :: Dx1, Dx2, Dy1, Dy2, Dz1, Dz2
    character (len=*), parameter :: NameSub=NameMod//'::trilinear_scalar'
    !--------------------------------------------------------------------------

    if ( present(iCell_D) .and. present(Dist_D) ) then

       ! Calculate the remaining cell indices and interpolation weights
       
       i1 = iCell_D(1)
       i2 = i1 + 1
       Dx1 = Dist_D(1)
       Dx2 = 1.0-Dx1

       j1 = iCell_D(2)
       j2 = j1 + 1
       Dy1 = Dist_D(2)
       Dy2 = 1.0 - Dy1

       k1 = iCell_D(3)
       k2 = k1 + 1
       Dz1 = Dist_D(3)
       Dz2 = 1.0 - Dz1

    else

       ! Locate the point Xyz_D on the grid
       call find_cell(iMin, iMax, Xyz_D(1), i1, Dx1, x_I, DoExtrapolate, &
            "Called for coord1 from "//NameSub)
       
       call find_cell(jMin, jMax, Xyz_D(2), j1, Dy1, y_I, DoExtrapolate, &
            "Called for coord2 from "//NameSub)

       call find_cell(kMin, kMax, Xyz_D(3), k1, Dz1, z_I, DoExtrapolate, &
            "Called for coord3 from "//NameSub)
    
       ! Calculate the remaining cell indices and interpolation weights
       i2 = i1 + 1; Dx2 = 1.0 - Dx1
       j2 = j1 + 1; Dy2 = 1.0 - Dy1
       k2 = k1 + 1; Dz2 = 1.0 - Dz1

    end if

    !Perform interpolation (or extrapolation)
    trilinear_scalar = Dz2*( Dy2*( Dx2*a_III(i1,j1,k1)   &
         +                         Dx1*a_III(i2,j1,k1))  &
         +                   Dy1*( Dx2*a_III(i1,j2,k1)   &
         +                         Dx1*a_III(i2,j2,k1))) &
         +             Dz1*( Dy2*( Dx2*a_III(i1,j1,k2)   &
         +                         Dx1*a_III(i2,j1,k2))  &
         +                   Dy1*( Dx2*a_III(i1,j2,k2)   &
         +                         Dx1*a_III(i2,j2,k2)))

  end function trilinear_scalar

  !===========================================================================

  function trilinear_vector( &
       a_VIII, nVar, iMin, iMax, jMin, jMax, kMin, kMax, Xyz_D, &
       x_I, y_I, z_I, DoExtrapolate, iCell_D, Dist_D)

    ! Calculate trilinear interpolation of a_III at position Xyz_D

    integer, intent(in) :: nVar, iMin, iMax, jMin, jMax, kMin, kMax
    real, intent(in)    :: a_VIII(nVar, iMin:iMax, jMin:jMax, kMin:kMax)
    real, intent(in)    :: Xyz_D(3)

    real,    intent(in), optional :: x_I(iMin:iMax)
    real,    intent(in), optional :: y_I(jMin:jMax)
    real,    intent(in), optional :: z_I(kMin:kMax)
    logical, intent(in), optional :: DoExtrapolate

    integer,    intent(in), optional :: iCell_D(3)
    real,    intent(in), optional :: Dist_D(3)

    ! return value
    real :: trilinear_vector(nVar)

    integer :: i1, i2, j1, j2, k1, k2
    real    :: Dx1, Dx2, Dy1, Dy2, Dz1, Dz2
    character (len=*), parameter :: NameSub=NameMod//'::trilinear_vector'
    !--------------------------------------------------------------------------

    if ( present(iCell_D) .and. present(Dist_D) ) then

       ! Calculate the remaining cell indices and interpolation weights

       i1 = iCell_D(1)
       i2 = i1 + 1
       Dx1 = Dist_D(1)
       Dx2 = 1.0 - Dx1

       j1 = iCell_D(2)
       j2 = j1 + 1
       Dy1 = Dist_D(2)
       Dy2 = 1.0 - Dy1

       k1 = iCell_D(3)
       k2 = k1 + 1
       Dz1 = Dist_D(3)
       Dz2 = 1.0 - Dz1

    else

       ! Locate the point Xyz_D on the grid
       call find_cell(iMin, iMax, Xyz_D(1), i1, Dx1, x_I, DoExtrapolate, &
            "Called for coord1 from "//NameSub)
       
       call find_cell(jMin, jMax, Xyz_D(2), j1, Dy1, y_I, DoExtrapolate, &
            "Called for coord2 from "//NameSub)

       call find_cell(kMin, kMax, Xyz_D(3), k1, Dz1, z_I, DoExtrapolate, &
            "Called for coord3 from "//NameSub)
    
       ! Calculate the remaining cell indices and interpolation weights
       i2 = i1 + 1; Dx2 = 1.0 - Dx1
       j2 = j1 + 1; Dy2 = 1.0 - Dy1
       k2 = k1 + 1; Dz2 = 1.0 - Dz1

    end if

    trilinear_vector = Dz2*(Dy2*(Dx2*a_VIII(:,i1,j1,k1)   &
         +                       Dx1*a_VIII(:,i2,j1,k1))  &
         +                  Dy1*(Dx2*a_VIII(:,i1,j2,k1)   &
         +                       Dx1*a_VIII(:,i2,j2,k1))) &
         +             Dz1*(Dy2*(Dx2*a_VIII(:,i1,j1,k2)   &
         +                       Dx1*a_VIII(:,i2,j1,k2))  &
         +                  Dy1*(Dx2*a_VIII(:,i1,j2,k2)   &
         +                       Dx1*a_VIII(:,i2,j2,k2)))

  end function trilinear_vector

  !===========================================================================
  subroutine find_cell(MinCoord, MaxCoord, Coord, iCoord, dCoord, &
       Coord_I, DoExtrapolate, StringError, IsInside)

    ! Find cell index and distance from cell for either 
    ! - a uniform grid with normalized coordinate (Coord_I is NOT present)
    ! - a nonuniform grid with monotone coordinates (Coord_I is present)
    !
    ! For sake of easy usage the returned coordinate index iCoord always
    ! satisfies MinCoord <= iCoord < MaxCoord.
    !
    ! If the coordinate is out of bounds, and DoExtrapolate is not present,
    ! the code stops with an error message. If DoExtrapolate is present and
    ! false, dCoord is modified to 0 or 1 so that the last grid cell is used.
    ! If DoExtrapolate is true, iCoord and dCoord are set 
    ! corresponding to linear extrapolation from the last two grid values.
    !
    ! For interpolation the normalized distance dCoord measured from 
    ! coordinate iCoord satisfies 0.0 <= dCoord <= 1.0
    ! but for extrapolation dCoord < 0.0 or > 1.0 is also possible.
    !---
    ! FOR THE UNIFORM CASE the normalized coordinate Coord should be equal to
    ! the index at the cell centers, therefore:
    !
    ! iCoord = max(MinCoord, min(MaxCoord-1, floor(Coord)))
    ! dCoord = Coord - iCoord
    !
    ! The optional IsInside = MinCoord <= Coord <= MaxCoord
    !----
    ! IN THE NON-UNIFORM CASE the cell iCoord that is left to coordinate Coord
    ! is found with a binary search in the Coord_I coordinates.
    !
    ! The normalized distance is set to
    !    dCoord = (Coord-Coord_I(iCoord))/(Coord_I(iCoord+1)-Coord_I(iCoord))
    !
    ! The optional IsInside = Coord_I(1) <= Coord <= Coord_I(nCoord).
    !----
    ! Example for linear interpolation on a 1D uniform grid of nX cells,
    ! DeltaX cell size and the first cell center is at DeltaX/2:
    !
    !   call find_cell(1, nX, x/DeltaX+0.5, iX, d)
    !   State_V = (1.0 - d)*State_VC(:,iX) + d*State_VC(:,iX+1)
    !
    ! Example for linear interpolation on a 1D non-uniform grid 
    ! with 2 ghost cells:
    !
    !   call find_cell(-1, nI+2, x, iX, d, x_G)
    !   State_V = (1.0 - d)*State_VG(:,iX) + d*State_VG(:,iX+1)

    integer, intent(in)           :: MinCoord, MaxCoord
    real,    intent(in)           :: Coord
    integer, intent(out)          :: iCoord
    real,    intent(out), optional:: dCoord
    real,    intent(in),  optional:: Coord_I(MinCoord:MaxCoord)
    logical, intent(in),  optional:: DoExtrapolate
    character(len=*),     optional:: StringError
    logical, intent(out), optional:: IsInside

    integer:: i, Di

    character(len=*), parameter:: NameSub="ModInterpolate::find_cell"
    !------------------------------------------------------------------------

    if(.not.present(Coord_I))then
       ! Uniform grid case with normalized coordinate

       iCoord = min(MaxCoord-1, max(MinCoord, floor(Coord)))
       dCoord = Coord - iCoord

       if(Coord < MinCoord)then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) NameSub, ': ', StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord=', Coord
             call CON_stop(NameSub//': normalized coordinate is to small!')
          elseif(.not.DoExtrapolate)then
             ! Use lefttmost cell (first order accurate) 
             dCoord = 0.0
          end if
          if(present(IsInside)) IsInside = .false.
       elseif(Coord > MaxCoord)then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord=', Coord
             call CON_stop(NameSub//': normalized coordinate is too large!')
          elseif(.not.DoExtrapolate)then
             ! Use rightmost cell (first order accurate) 
             dCoord = 1.0
          endif
          if(present(IsInside)) IsInside = .false.
       else
          if(present(IsInside)) IsInside = .true.
       end if

    elseif(Coord_I(MinCoord) < Coord_I(MaxCoord))then

       ! Monotone increasing coordinates

       if(Coord < Coord_I(MinCoord))then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord, CoordMin, CoordMax=', &
                  Coord, Coord_I(MinCoord), Coord_I(MaxCoord)
             call CON_stop(NameSub//': coordinate is too small!')
          elseif(DoExtrapolate)then
             iCoord = MinCoord
             dCoord = (Coord - Coord_I(iCoord)) &
                  /   (Coord_I(iCoord+1) - Coord_I(iCoord))
          else
             iCoord = MinCoord
             dCoord = 0.0
          end if
          if(present(IsInside)) IsInside = .false.
          RETURN
       end if

       if(Coord > Coord_I(MaxCoord))then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord, CoordMin, CoordMax=', &
                  Coord, Coord_I(MinCoord), Coord_I(MaxCoord)
             call CON_stop(NameSub//': coordinate is too large!')
          elseif(DoExtrapolate)then
             iCoord = MaxCoord - 1
             dCoord = (Coord - Coord_I(iCoord))  &
                  /   (Coord_I(iCoord+1) - Coord_I(iCoord))
          else
             iCoord = MaxCoord - 1
             dCoord = 1.0
          end if
          if(present(IsInside)) IsInside = .false.
          RETURN
       end if

       if(present(IsInside)) IsInside = .true.

       ! binary search
       i  = (MinCoord + MaxCoord)/2
       Di = (MaxCoord - MinCoord)/2
       do
          Di = (Di + 1)/2
          if(Coord < Coord_I(i)) then
             i = max(MinCoord, i - Di)
          elseif(Coord > Coord_I(i+1))then
             i = min(MaxCoord-1, i + Di)
          else
             EXIT
          end if
       end do
       iCoord = i
       if(Coord_I(iCoord+1) == Coord_I(iCoord))then
          dCoord = 0.0
       else
          dCoord = (Coord             - Coord_I(iCoord)) &
               /   (Coord_I(iCoord+1) - Coord_I(iCoord))
       end if
    else

       ! Monotone decreasing coordinates

       if(Coord < Coord_I(MaxCoord))then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord, CoordMin, CoordMax=', &
                  Coord, Coord_I(MaxCoord), Coord_I(MinCoord)
             call CON_stop(NameSub//': coordinate is too small!')
          elseif(DoExtrapolate)then
             iCoord = MaxCoord - 1
             dCoord = (Coord_I(iCoord) - Coord) &
                  /   (Coord_I(iCoord) - Coord_I(iCoord+1))
          else
             iCoord = MaxCoord - 1
             dCoord = 1.0
          end if
          if(present(IsInside)) IsInside = .false.
          RETURN
       end if

       if(Coord > Coord_I(MinCoord))then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord, CoordMin, CoordMax=', &
                  Coord, Coord_I(MaxCoord), Coord_I(MinCoord)
             call CON_stop(NameSub//': coordinate is too large!')
          elseif(DoExtrapolate)then
             iCoord = MinCoord
             dCoord = (Coord_I(iCoord) - Coord)  &
                  /   (Coord_I(iCoord) - Coord_I(iCoord+1))
          else
             iCoord = MinCoord
             dCoord = 0.0
          end if
          if(present(IsInside)) IsInside = .false.
          RETURN
       end if

       if(present(IsInside)) IsInside = .true.

       ! binary search
       i  = (MinCoord + MaxCoord)/2
       Di = (MaxCoord - MinCoord)/2
       do
          Di = (Di + 1)/2
          if(Coord > Coord_I(i)) then
             i = max(MinCoord, i - Di)
          elseif(Coord < Coord_I(i+1))then
             i = min(MaxCoord-1, i + Di)
          else
             EXIT
          end if
       end do
       iCoord = i
       if(Coord_I(iCoord+1) == Coord_I(iCoord))then
          dCoord = 0.0
       else
          dCoord = (Coord_I(iCoord) - Coord  ) &
               /   (Coord_I(iCoord) - Coord_I(iCoord+1))
       end if
    end if

  end subroutine find_cell
  !===========================================================================
  subroutine fit_parabola(x_I, y_I, &
       xExtremumOut, yExtremumOut, Weight2Out_I, Weight3Out_I)

    ! Given 3 discrete points at x_D and 3 function values y_D
    ! with the middle point being the discrete extrem value,
    ! find the extremum value of the parabola going through the points,
    ! and the 3rd order interpolation weights to interpolate to this point

    real, intent(in)           :: x_I(3)         ! coordinates
    real, intent(in)           :: y_I(3)         ! values
    real, intent(out), optional:: xExtremumOut   ! coordinate of extremum
    real, intent(out), optional:: yExtremumOut   ! value of extremum
    real, intent(out), optional:: Weight2Out_I(3)! weights for 2nd order interpolation
    real, intent(out), optional:: Weight3Out_I(3)! weights for 3rd order interpolation

    real:: xE, yE          ! coordinates of extremum
    real:: x1, y1, x3, y3  ! shifted coordinates of points 1 and 3
    real:: s1, s3          ! slopes of 1-2 and 2-3 segments

    real:: Ratio, Area2

    character(len=*), parameter:: NameSub = 'fit_parabola'
    !------------------------------------------------------------------------

    ! Shift coordinates so that x2 = 0
    x1 = x_I(1) - x_I(2)
    x3 = x_I(3) - x_I(2)

    ! Shift values so that y2 = 0
    y1 = y_I(1) - y_I(2)
    y3 = y_I(3) - y_I(2)

    if(x1 == 0.0 .or. x3 == 0.0)then
       write(*,*) NameSub,': x_I=', x_I,' y_I=', y_I
       call CON_stop(NameSub//' error in coordinates')
    end if

    ! Calculate slopes
    s1 = y1/x1
    s3 = y3/x3

    if(s1*s3 > 0.0)then
       write(*,*) NameSub,': x_I=', x_I,' y_I=', y_I
       call CON_stop(NameSub//' error: midpoint is not an extremum')
    end if

    ! Find the position where the line connecting 
    ! the (x1/2, s1) and (x3/2, s3) points intersects the X axis.
    ! This is where the slope of the parabola is zero

    xE = 0.5*x1 + 0.5*(x3 - x1)*s1/(s1 - s3)

    if(present(xExtremumOut)) xExtremumOut = xE + x_I(2)

    if(present(Weight2Out_I))then
       ! Use the two points surrounding xE for linear interpolation
       if(xE > 0.0) then
          Weight2Out_I(1) = 0.0
          Weight2Out_I(3) = xE/x3
          Weight2Out_I(2) = 1.0 - Weight2Out_I(3)
       else
          Weight2Out_I(3) = 0.0
          Weight2Out_I(1) = xE/x1
          Weight2Out_I(2) = 1.0 - Weight2Out_I(1)
       end if
    end if

    if(present(yExtremumOut) .or. present(Weight3Out_I))then
       ! Find the value of the parabola y = a*(x-xE)**2 + yE at the extremum
       ! We can use any 2 of the points to solve for yE.
       if(xE > 0.0)then
          Ratio = xE**2/(xE - x1)**2
          yE = Ratio*y1/(Ratio - 1.0)
       else
          Ratio = xE**2/(x3 - xE)**2
          yE = Ratio*y3/(Ratio - 1.0)
       end if
       if(present(yExtremumOut)) yExtremumOut = yE + y_I(2)

       if(present(Weight3Out_I))then

          ! Calculate 3rd order interpolation weights from the 3 points 
          ! to the location of the extremum.
          ! We use the fact that the parabola is an exact solution.

          ! Twice the area of the triangle with sign  (x1,y1) x (x3,y3)
          Area2 = x1*y3 - y1*x3

          ! For points 1 and 3 the weight is the fraction of the triangle
          ! Area(2,3,E)/Area(1,2,3) and Area(1,2,E)/Area(1,2,3)
          Weight3Out_I(1) =  (xE*y3 - yE*x3)/Area2
          Weight3Out_I(3) = -(xE*y1 - yE*x1)/Area2

          ! For point 2 we use that the sum of weights must be 1
          Weight3Out_I(2) = 1.0 - Weight3Out_I(1) - Weight3Out_I(3)
       end if
    end if

  end subroutine fit_parabola
  !===========================================================================
  subroutine test_interpolation

    real :: a_I(0:2) = (/ 10., 20., 30. /)

    real :: a_II(2,3) = reshape((/ 1., 20., 3., 40., 5., 60. /), (/2, 3/))

    real :: a_VII(2,2,3) = reshape( &
         (/1., 10., 20., 200., 3., 30., 40., 400., 5., 50., 60., 600./), &
         (/2, 2, 3/))

    real :: a_III(2,2,0:2) = reshape((/ &
         1., 20., 3., 40., &
         100., 2000., 300., 4000., &
         10000., 200000., 30000., 400000. /), (/2, 2, 3/))

    real :: a_VIII(2,2,2,0:2) = reshape((/ &
         1., -10., 20., -200., 3., -30., 40., -400., &
         100., -1000., 2000., -20000., 300., -3000., 4000., -40000.,  &
         1e4, -1e5, 2e5, -2e6, 3e4, -3e5, 4e5, -4e6 /), (/2, 2, 2, 3/))

    real :: x12_I(1:2) = (/ 1., 2./)
    real :: x13_I(1:3) = (/ 1., 2., 4./)
    real :: x02_I(0:2) = (/ 1., 2., 4./)

    integer, parameter:: MinCoord = 1, MaxCoord = 8
    real   :: Coord_I(MinCoord:MaxCoord) = &
         (/ 1.0, 2.0, 4.0, 8.0, 16.0, 17.0, 24.0, 25.0 /)
    integer:: nCoord, iSign
    real   :: Coord, dCoord, CoordMin, CoordMax
    integer:: i, iCoord
    logical:: IsInside

    real :: Result, GoodResult, Result_V(2), GoodResult_V(2)

    ! Variables for fit_parabola test
    real:: x_I(3), y_I(3), xMin, yMin, xExtremum, yExtremum
    real:: Weight2_I(3), Weight3_I(3)

    character(len=*), parameter:: NameSub=NameMod//"::test_interpolation"
    !----------------------------------------------------------------------
    ! Change sign of coordinates to test for increasing and decreasing orders
    do iSign = 1, -1, -2
       if(iSign == 1)then
          write(*,'(a)')'Testing find_cell for increasing coordinates'
       else
          write(*,'(a)')'Testing find_cell for decreasing coordinates'
       end if

       ! Change number of coordinates to test binary search
       do nCoord = MaxCoord/2, MaxCoord

          ! Search for all integer coordinates
          ! starting below and finishing above the coordinate range

          CoordMin = min(Coord_I(MinCoord), Coord_I(nCoord))
          CoordMax = max(Coord_I(MinCoord), Coord_I(nCoord))

          do i = ceiling(CoordMin) - 1, floor(CoordMax) + 1
             Coord = i
             call find_cell(MinCoord, nCoord, Coord, &
                  iCoord, dCoord, Coord_I, .false., &
                  'Called from '//NameSub, IsInside)

             if(iSign*Coord < iSign*Coord_I(MinCoord))then
                if(IsInside) write(*,*) &
                     'Test failed for nCoord, Coord=', nCoord, Coord, &
                     ', IsInside=T, should be false'
                if(iCoord /= MinCoord) write(*,*)&
                     'Test failed for nCoord, Coord=', nCoord, Coord, &
                     ', iCoord=', iCoord, ' should be ', MinCoord
                if(dCoord /= 0.0) write(*,*) &
                     'Test failed for nCoord, Coord=', nCoord, Coord, &
                     ', dCoord=', dCoord, ' should be 0.0'
                CYCLE
             end if
             if(iSign*Coord > iSign*Coord_I(nCoord))then
                if(IsInside) write(*,*) &
                     'Test failed for nCoord, Coord=', nCoord, Coord, &
                     ', IsInside=T, should be false'
                if(iCoord /= nCoord - 1) write(*,*) &
                     'Test failed for nCoord, Coord=', nCoord, Coord, &
                     ', iCoord=', iCoord, ' should be ', nCoord - 1
                if(dCoord /= 1.0) write(*,*) &
                     'Test failed for nCoord, Coord=', nCoord, Coord, &
                     ', dCoord=', dCoord, ' should be 1.0'
                CYCLE
             end if
             if(.not.IsInside) write(*,*) &
                  'Test failed for nCoord, Coord=', nCoord, Coord, &
                  ', IsInside=F, should be true'

             if(iCoord < MinCoord .or. iCoord > nCoord-1) then
                write(*,*) &
                     'Test failed for nCoord, Coord=', nCoord, Coord, &
                     ', iCoord=', iCoord, ' should be < ', MinCoord, &
                     ' and > ', nCoord - 1
                CYCLE
             end if

             if(iSign*Coord_I(iCoord) > iSign*Coord) write(*,*) &
                  'Test failed for nCoord, Coord=', nCoord, Coord, &
                  ', iSign*Coord_I(iCoord)=', iSign*Coord_I(iCoord), &
                  ' should be <= iSign*Coord'

             if(iSign*Coord_I(iCoord+1) < iSign*Coord) write(*,*)       &
                  'Test failed for nCoord, Coord=', nCoord, Coord, &
                  ', iSign*Coord_I(iCoord+1)=', iSign*Coord_I(iCoord+1), &
                  ' should be >= iSign*Coord' 
             if(abs(Coord_I(iCoord) &
                  + dCoord*(Coord_I(iCoord+1) - Coord_I(iCoord)) &
                  - Coord) > 1e-6) write(*,*) &
                  'Test failed for nCoord, Coord=', nCoord, Coord, &
                  ', Coord_I(iCoord:iCoord+1)=', Coord_I(iCoord:iCoord+1), &
                  ', but incorrect dCoord = ', dCoord
          end do
       end do
       ! Change signs of coordinates to test decreasing order
       Coord_I = -Coord_I
    end do

    !Test for normal conditions.
    write(*,'(a)')'Testing function linear for uniform grid'
    Result = linear(a_I, 0, 2, 1.1)
    GoodResult = 21.0
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing function linear for non-uniform grid'
    Result = linear(a_I, 0, 2, 2.2, x02_I)
    GoodResult = 21.0
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing function bilinear for uniform grid'
    Result = bilinear(a_II, 1, 2, 1, 3, (/1.1, 2.2/))
    GoodResult = 7.46
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing function bilinear for non-uniform grid'
    Result = bilinear(a_II, 1, 2, 1, 3, (/1.1, 2.2/), x12_I, x13_I)
    GoodResult = 7.08
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing function trilinear for uniform grid'
    Result = trilinear(a_III, 1, 2, 1, 2, 0, 2, (/1.1, 1.2, 1.3/))
    GoodResult = 11236.2
    if(abs(Result - GoodResult) > 1.e-2) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing function trilinear for nonuniform grid'
    Result = trilinear(a_III, 1, 2, 1, 2, 0, 2, (/1.1, 1.2, 1.3/), &
         x12_I, x12_I, x02_I)
    GoodResult = 112.362
    if(abs(Result - GoodResult) > 1.e-2) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    !Test out-of-bounds, no extrapolation
    write(*,'(a)')'Testing bilinear out-of-bounds: +X for uniform grid'
    Result = bilinear(a_II, 1, 2, 1, 3, (/3.,1./), DoExtrapolate=.false.)
    GoodResult = 20.0
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing bilinear out-of-bounds: +X for nonuniform grid'
    Result = bilinear(a_II, 1, 2, 1, 3, (/3.,1./), x12_I, x13_I, &
         DoExtrapolate=.false.)
    GoodResult = 20.0
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing bilinear out-of-bounds: -X for uniform grid'
    Result = bilinear(a_II, 1, 2, 1, 3, (/-3.,2./), DoExtrapolate=.false.)
    GoodResult = 3.0
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing bilinear out-of-bounds: -X for nonuniform grid'
    Result = bilinear(a_II, 1, 2, 1, 3, (/-3.,2./), x12_I, x13_I, &
         DoExtrapolate=.false.)
    GoodResult = 3.0
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing bilinear out-of-bounds: +Y for uniform grid'
    Result = bilinear(a_II, 1, 2, 1, 3, (/1.,6./), DoExtrapolate=.false.)
    GoodResult = 5.0
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing bilinear out-of-bounds: +Y for nonuniform grid'
    Result = bilinear(a_II, 1, 2, 1, 3, (/1.,6./), x12_I, x13_I, &
         DoExtrapolate=.false.)
    GoodResult = 5.0
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing bilinear out-of-bounds: -Y for uniform grid'
    Result = bilinear(a_II, 1, 2, 1, 3, (/2.,-3./), DoExtrapolate=.false.)
    GoodResult = 20.0
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing bilinear out-of-bounds: -Y for nonuniform grid'
    Result = bilinear(a_II, 1, 2, 1, 3, (/2.,-3./), x12_I, x13_I, &
         DoExtrapolate=.false.)
    GoodResult = 20.0
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing trilinear out-of-bounds: +Z for uniform grid'
    Result = trilinear(a_III, 1, 2, 1, 2, 0, 2, (/1., 1., 2.4/), &
         DoExtrapolate=.false.)
    GoodResult = 10000.0
    if(abs(Result - GoodResult) > 1.e-2) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing trilinear out-of-bounds: +Z for nonuniform grid'
    Result = trilinear(a_III, 1, 2, 1, 2, 0, 2, (/1., 1., 4.1/), &
         x12_I, x12_I, x02_I, DoExtrapolate=.false.)
    GoodResult = 10000.0
    if(abs(Result - GoodResult) > 1.e-6) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing trilinear out-of-bounds: -Z for uniform grid'
    Result = trilinear(a_III, 1, 2, 1, 2, 0, 2, (/1., 1., -0.4/), &
         DoExtrapolate=.false.)
    GoodResult = 1.0
    if(abs(Result - GoodResult) > 1.e-6) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing trilinear out-of-bounds: -Z for nonuniform grid'
    Result = trilinear(a_III, 1, 2, 1, 2, 0, 2, (/1., 1., 0.1/), &
         x12_I, x12_I, x02_I, DoExtrapolate=.false.)
    GoodResult = 1.0
    if(abs(Result - GoodResult) > 1.e-2) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    !Test extrapolation
    write(*,'(a)')'Testing bilinear extrapolation: +X uniform'
    Result = bilinear(a_II, 1, 2, 1, 3, (/2.5,1./), DoExtrapolate=.true.)
    GoodResult = 29.5
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult
    
    write(*,'(a)')'Testing bilinear extrapolation: +X nonuniform'
    Result = bilinear(a_II, 1, 2, 1, 3, (/2.5,1./), x12_I, x13_I, &
         DoExtrapolate=.true.)
    GoodResult = 29.5
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult
    
    write(*,'(a)')'Testing bilinear extrapolation: -X uniform'
    Result = bilinear(a_II, 1, 2, 1, 3, (/.5,1.5/), DoExtrapolate=.true.)
    GoodResult = -12.0
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing bilinear extrapolation: -X nonuniform'
    Result = bilinear(a_II, 1, 2, 1, 3, (/.5,1.5/), x12_I, x13_I, &
         DoExtrapolate=.true.)
    GoodResult = -12.0
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing trilinear extrapolation: +Z uniform'
    Result = trilinear(a_III, 1, 2, 1, 2, 0, 2, (/1.3, 1.9, 2.60/), &
         DoExtrapolate=.true.)
    GoodResult = 212958.38
    if(abs(Result - GoodResult) > 1.0) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing trilinear extrapolation: +Z nonuniform'
    Result = trilinear(a_III, 1, 2, 1, 2, 0, 2, (/1.3, 1.9, 5.2/), &
         x12_I, x12_I, x02_I, DoExtrapolate=.true.)
    GoodResult = 212958.38
    if(abs(Result - GoodResult) > 1.0) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing function bilinear_vector'
    Result_V = bilinear(a_VII, 2, 1, 2, 1, 3, (/1.1, 2.2/))
    GoodResult_V = (/7.46, 74.6/)
    if(any(abs(Result_V - GoodResult_V) > 1.e-5)) &
         write(*,*) 'Test failed: Result=',Result_V,&
         ' differs from ',GoodResult_V

    write(*,'(a)')'Testing function trilinear_vector'
    Result_V = trilinear(a_VIII, 2, 1, 2, 1, 2, 0, 2, (/1.1, 1.2, 1.3/))
    GoodResult_V = (/ 11236.2, -112362.0 /)
    if(any(abs(Result_V - GoodResult_V) > 1.e-2)) write(*,*) &
         'Test failed: Result=', Result_V, ' differs from ', GoodResult_V

    write(*,'(a)')'Testing fit_parabola'
    x_I = (/ 3.1, 4.0, 5.0 /)
    xMin = 4.2; yMin = 1.5
    y_I = 0.1*(x_I - xMin)**2 + yMin
    call fit_parabola(x_I, y_I, xExtremum, yExtremum, Weight2_I, Weight3_I)

    ! write(*,*)'x_I, xE=', x_I, xExtremum
    ! write(*,*)'y_I, yE=', y_I, yExtremum
    ! write(*,*)'Weight2_I=', Weight2_I
    ! write(*,*)'Weight3_I=', Weight3_I

    if(abs(xExtremum - xMin) > 1e-6) write(*,*) &
         'Test failed: xExtremum=', xExtremum, ' differs from ', xMin

    if(abs(yExtremum - yMin) > 1e-6) write(*,*) &
         'Test failed: yExtremum=', yExtremum, ' differs from ', yMin

    if(abs(sum(Weight2_I) - 1.0) > 1e-6) write(*,*) &
         'Test failed: sum of Weight2_I=', Weight2_I, ' is not 1'

    if(abs(sum(Weight2_I*x_I) - xMin) > 1e-6) write(*,*) &
         'Test failed: Weight2_I=', Weight2_I, ' sum(Weight2_I*x_I)=', &
         sum(Weight2_I*x_I), ' differs from ', xMin

    if(abs(sum(Weight3_I) - 1.0) > 1e-6) write(*,*) &
         'Test failed: sum of Weight3_I=', Weight3_I, ' is not 1'

    if(abs(sum(Weight3_I*x_I) - xMin) > 1e-6) write(*,*) &
         'Test failed: Weight3_I=', Weight3_I, ' sum(Weight3_I*x_I)=', &
         sum(Weight3_I*x_I), ' differs from ', xMin

    if(abs(sum(Weight3_I*y_I) - yMin) > 1e-6) write(*,*) &
         'Test failed: Weight3_I=', Weight3_I, ' sum(Weight3_I*y_I)=', &
         sum(Weight3_I*y_I), ' differs from ', yMin
    
  end subroutine test_interpolation

end module ModInterpolate
