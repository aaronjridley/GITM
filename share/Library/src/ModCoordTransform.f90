!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!BOP -------------------------------------------------------------------
!
!MODULE: ModCoordTransform - rotation matrices, spherical/Cartesian transforms
!
!DESCRIPTION:
!
! Contains general subroutines and functions to construct 3 by 3 rotation
! matrices around the principal axes, rotation matrix that transform between
! Cartesian and spherical vectors, and subroutines that convert between
! Cartesian and spherical coordinates and directions.
!
! All functions and subroutines have several versions so that they can
! be called with various arguments:
!
! A Cartesian position or vector can be given as
! \begin{itemize}
! \item 1 array with 3 elements in x, y, z order
! \item 3 scalars in x, y, z order
! \end{itemize}
!
! A spherical position or vector can be given as
! \begin{itemize}
! \item 1 array with 3 elements in r, $\theta$, $\phi$ order
! \item 3 scalars in r, $\theta$, $\phi$ order.
! \end{itemize}
!
! A direction can be given as
! \begin{itemize}
! \item 1 Cartesian vector with non-zero length
! \item 3 Cartesian components in x, y, z order
! \item 2 spherical angles in the $\theta$, $\phi$ order
! \item 4 scalars in the $\sin\theta, \cos\theta, \sin\phi, \cos\phi$ order
! \end{itemize}
!
! A rotational angle can be given as
! \begin{itemize}
! \item 1 scalar angle $\alpha$ in radians
! \item 2 scalars in the $\sin\alpha, \cos\alpha$ order.
! \end{itemize}
!
! As a convenient utility, the {\bf cross\_product} function is provided.
! It returns the cross product of the 2 Cartesian vectors as a 3 element array.
!INTERFACE:

module ModCoordTransform

  use ModNumConst, ONLY: cTwoPi, cUnit_DD

  implicit none

  save

  private ! except

  !PUBLIC MEMBER FUNCTIONS:

  public :: rot_matrix_x   ! rotation matrix around X axis (angle)
  public :: rot_matrix_y   ! rotation matrix around Y axis (angle)
  public :: rot_matrix_z   ! rotation matrix around Z axis (angle)
  public :: rot_xyz_sph    ! rotation matrix between Cartesian-spherical (dir)
  public :: xyz_to_sph     ! convert Cartesian into spherical coordinates 
  public :: sph_to_xyz     ! convert spherical into Cartesian coordinates
  public :: xyz_to_dir     ! convert Cartesian vector to spherical direction
  public :: dir_to_xyz     ! convert spher. direction to Cartesian unit vector
  public :: cross_product  ! return the cross product of two vectors
  public :: inverse_matrix ! return the inverse of a 3 by 3 matrix

  public :: show_rot_matrix      ! write out matrix elements in a nice format
  public :: atan2_check          ! compute atan2 even if both x and y are zero
  public :: test_coord_transform ! unit tester

  !REVISION HISTORY:
  ! 08Aug03 - Gabor Toth <gtoth@umich.edu> - initial prototype/prolog/code
  ! 29Jun06 - YingJuan - added inverse_matrix function
  !EOP ___________________________________________________________________

  integer, parameter :: x_=1, y_=2, z_=3

  interface rot_matrix_x
     module procedure rot_matrix_x1, rot_matrix_x2
  end interface

  interface rot_matrix_y
     module procedure rot_matrix_y1, rot_matrix_y2
  end interface

  interface rot_matrix_z
     module procedure rot_matrix_z1, rot_matrix_z2
  end interface

  interface rot_xyz_sph
     module procedure rot_xyz_sph1, rot_xyz_sph2, rot_xyz_sph3, rot_xyz_sph4
  end interface

  interface xyz_to_sph
     module procedure xyz_to_sph11, xyz_to_sph13, xyz_to_sph31, xyz_to_sph33
  end interface

  interface sph_to_xyz
     module procedure sph_to_xyz11, sph_to_xyz13, sph_to_xyz31, sph_to_xyz33
  end interface

  interface xyz_to_dir
     module procedure xyz_to_dir12, xyz_to_dir32, xyz_to_dir14, xyz_to_dir34
  end interface

  interface dir_to_xyz
     module procedure dir_to_xyz21, dir_to_xyz23, dir_to_xyz41, dir_to_xyz43
  end interface

  interface cross_product
     module procedure cross_product11, cross_product13, cross_product31,&
          cross_product33
  end interface

  character(len=*), parameter :: NameMod='ModCoordTransform'

contains

  !============================================================================
  subroutine xyz_to_sph11(Xyz_D,Sph_D)

    real, intent(in) :: Xyz_D(3)
    real, intent(out):: Sph_D(3)

    call xyz_to_sph(Xyz_D(1),Xyz_D(2),Xyz_D(3),Sph_D(1),Sph_D(2),Sph_D(3))

  end subroutine xyz_to_sph11

  !============================================================================
  subroutine xyz_to_sph13(Xyz_D,r,Theta,Phi)

    real, intent(in) :: Xyz_D(3)
    real, intent(out):: r,Theta,Phi

    call xyz_to_sph(Xyz_D(1),Xyz_D(2),Xyz_D(3),r,Theta,Phi)

  end subroutine xyz_to_sph13

  !============================================================================
  subroutine xyz_to_sph31(x,y,z,Sph_D)

    real, intent(in) :: x,y,z
    real, intent(out):: Sph_D(3)

    call xyz_to_sph(x,y,z,Sph_D(1),Sph_D(2),Sph_D(3))

  end subroutine xyz_to_sph31

  !BOP -------------------------------------------------------------------
  !IROUTINE: xyz_to_sph - convert Cartesian into spherical coordinates
  !
  !INTERFACE:

  subroutine xyz_to_sph33(x,y,z,r,Theta,Phi)

    !INPUT PARAMETERS:
    real, intent(in) :: x,y,z

    !OUTPUT PARAMETERS:
    real, intent(out):: r,Theta,Phi

    !DESCRIPTION:
    !
    ! This is the fundamental routine that converts Cartesian position into
    ! spherical position. All other xyz\_to\_sph varieties call this 
    ! subroutine. Here both vectors are given with 3 scalars 
    ! (hence the 33 suffix).

    !LOCAL VARIABLES:
    real :: d
    !EOP 
    !___________________________________________________________________
    !BOC
    d     = x**2 + y**2
    r     = sqrt(d + z**2)
    d     = sqrt(d)
    Theta = atan2_check(d,z)
    Phi   = atan2_check(y,x)
    !EOC
  end subroutine xyz_to_sph33

  !============================================================================

  subroutine sph_to_xyz11(Sph_D,Xyz_D)

    real, intent(in) :: Sph_D(3)
    real, intent(out):: Xyz_D(3)

    call sph_to_xyz(Sph_D(1),Sph_D(2),Sph_D(3),Xyz_D(1),Xyz_D(2),Xyz_D(3))

  end subroutine sph_to_xyz11

  !============================================================================

  subroutine sph_to_xyz31(r,Theta,Phi,Xyz_D)

    real, intent(in) :: r,Theta,Phi
    real, intent(out):: Xyz_D(3)

    call sph_to_xyz(r,Theta,Phi,Xyz_D(1),Xyz_D(2),Xyz_D(3))

  end subroutine sph_to_xyz31

  !============================================================================

  subroutine sph_to_xyz13(Sph_D,x,y,z)

    real, intent(in) :: Sph_D(3)
    real, intent(out):: x,y,z

    call sph_to_xyz(Sph_D(1),Sph_D(2),Sph_D(3),x,y,z)

  end subroutine sph_to_xyz13

  !BOP -------------------------------------------------------------------
  !IROUTINE: sph_to_xyz - convert spherical into Cartesian coordinates
  !
  !INTERFACE:

  subroutine sph_to_xyz33(r,Theta,Phi,x,y,z)

    !INPUT PARAMETERS:
    real, intent(in)  :: r,Theta,Phi

    !OUTPUT PARAMETERS:
    real, intent(out) :: x,y,z

    !DESCRIPTION:
    ! This is the fundamental routine that converts spherical position into
    ! Cartesian position. All other sph\_to\_xyz varieties call this
    ! subroutine. Here both vectors are given with 3 scalars
    ! (hence the 33 suffix).

    !LOCAL VARIABLES:
    real :: SinTheta
    !EOP
    !------------------------------------------------------------------------
    !BOC
    SinTheta = sin(Theta)
    x = r*SinTheta*cos(Phi)
    y = r*SinTheta*sin(Phi)
    z = r*cos(Theta)
    !EOC

  end subroutine sph_to_xyz33

  !============================================================================

  subroutine xyz_to_dir12(Xyz_D,Theta,Phi)

    real, intent(in) :: Xyz_D(3)
    real, intent(out):: Theta,Phi

    call xyz_to_dir(Xyz_D(1),Xyz_D(2),Xyz_D(3),Theta,Phi)

  end subroutine xyz_to_dir12

  !BOP -------------------------------------------------------------------
  !IROUTINE: xyz_to_dir - convert Cartesian vector into direction (angles)
  !
  !INTERFACE:
  subroutine xyz_to_dir32(x,y,z,Theta,Phi)

    !INPUT PARAMETERS:
    real, intent(in) :: x,y,z

    !OUTPUT PARAMETERS:
    real, intent(out):: Theta,Phi

    !DESCRIPTION:
    !
    ! This is the fundamental routine that converts Cartesian vector into
    ! spherical direction given with angles. 
    ! The Cartesian vector can have any positive length.
    !EOP
    !-----------------------------------------------------------------------
    !BOC

    Theta = atan2_check(sqrt(x**2 + y**2),z)
    Phi   = atan2_check(y,x)

    !EOC
  end subroutine xyz_to_dir32

  !============================================================================

  subroutine xyz_to_dir14(Xyz_D,SinTheta,CosTheta,SinPhi,CosPhi)

    real, intent(in) :: Xyz_D(3)
    real, intent(out):: SinTheta,CosTheta,SinPhi,CosPhi

    call xyz_to_dir(Xyz_D(1),Xyz_D(2),Xyz_D(3),SinTheta,CosTheta,SinPhi,CosPhi)

  end subroutine xyz_to_dir14

  !============================================================================

  !BOP -------------------------------------------------------------------
  !IROUTINE: xyz_to_dir - convert Cartesian vector into direction (trigonometric)
  !
  !INTERFACE:
  subroutine xyz_to_dir34(x,y,z,SinTheta,CosTheta,SinPhi,CosPhi)

    !INPUT PARAMETERS:
    real, intent(in) :: x,y,z

    !OUTPUT PARAMETERS:
    real, intent(out):: SinTheta,CosTheta,SinPhi,CosPhi

    !DESCRIPTION:
    !
    ! This is the fundamental routine that converts Cartesian vector into
    ! spherical direction given with trigonometric functions.
    ! This version is faster (in principle) than the one using angles.
    ! On the other hand here one needs to check the special case 
    ! when the vector is parallel with the Z axis.
    ! The Cartesian vector can have any positive length.

    !LOCAL VARIABLES:
    real :: r,d
    !EOP
    !-----------------------------------------------------------------------
    !BOC
    d = x**2+y**2
    r = sqrt(d+z**2)
    d = sqrt(d)
    SinTheta = d/r
    CosTheta = z/r
    if(d > 0)then
       SinPhi   = y/d
       CosPhi   = x/d
    else
       SinPhi   = 0
       CosPhi   = 0
    end if
    !EOC

  end subroutine xyz_to_dir34

  !============================================================================

  subroutine dir_to_xyz21(Theta,Phi,Xyz_D)

    real, intent(in) :: Theta,Phi
    real, intent(out):: Xyz_D(3)

    call dir_to_xyz(sin(Theta),cos(Theta),sin(Phi),cos(Phi),&
         Xyz_D(1),Xyz_D(2),Xyz_D(3))

  end subroutine dir_to_xyz21

  !============================================================================

  subroutine dir_to_xyz23(Theta,Phi,x,y,z)

    real, intent(in) :: Theta,Phi
    real, intent(out):: x,y,z
    
    call dir_to_xyz(sin(Theta),cos(Theta),sin(Phi),cos(Phi),x,y,z)

  end subroutine dir_to_xyz23

  !============================================================================

  subroutine dir_to_xyz41(SinTheta,CosTheta,SinPhi,CosPhi,Xyz_D)

    real, intent(in) :: SinTheta,CosTheta,SinPhi,CosPhi
    real, intent(out):: Xyz_D(3)

    call dir_to_xyz(SinTheta,CosTheta,SinPhi,CosPhi,Xyz_D(1),Xyz_D(2),Xyz_D(3))

  end subroutine dir_to_xyz41

  !BOP -------------------------------------------------------------------
  !IROUTINE: dir_to_xyz - convert spher. direction into Cartesian unit vector
  !
  !INTERFACE:
  subroutine dir_to_xyz43(SinTheta,CosTheta,SinPhi,CosPhi,x,y,z)

    !INPUT PARAMETERS:
    real, intent(in) :: SinTheta,CosTheta,SinPhi,CosPhi

    !OUTPUT PARAMETERS:
    real, intent(out):: x,y,z

    !DESCRIPTION:
    !
    ! This is the fundamental routine that converts a spherical direction 
    ! into a Cartesian unit vector. The spherical direction can be
    ! given with 4 trigonometric functions or 2 angles.
    ! The Cartesian unit vector can be returned in 1 array or 3 scalars.

    !EOP
    !-----------------------------------------------------------------------
    !BOC
    x = SinTheta*CosPhi
    y = SinTheta*SinPhi
    z = CosTheta
    !EOC

  end subroutine dir_to_xyz43

  !============================================================================

  function rot_matrix_x1(Angle) result(Rot_DD)

    real, intent(in) :: Angle
    real :: Rot_DD(3,3)

    rot_DD = rot_matrix_x(sin(Angle),cos(Angle))

  end function rot_matrix_x1

  !BOP -------------------------------------------------------------------
  !IROUTINE: rot_matrix_x - rotation matrix around X axis
  !
  !INTERFACE:
  function rot_matrix_x2(SinAngle, CosAngle) result(Rot_DD)

    !INPUT PARAMETERS:
    real, intent(in) :: SinAngle, CosAngle

    !RETURN VALUE:
    real :: Rot_DD(3,3)

    !DESCRIPTION:
    !
    ! This is the fundamental routine that calculates the rotation
    ! matrix around the X axis. Here the rotation angle is given with 2 
    ! trigonometric functions. Alternatively it can be given with 1 angle
    ! in radians. The rot\_matrix\_y and rot\_matrix\_z functions are
    ! fully analogous.
    !EOP
    !----------------------------------------------------------------------
    !BOC
    Rot_DD        =  0
    Rot_DD(y_,y_) =  CosAngle
    Rot_DD(y_,z_) = -SinAngle
    Rot_DD(z_,y_) =  SinAngle
    Rot_DD(z_,z_) =  CosAngle
    Rot_DD(x_,x_) =  1
    !EOC
  end function rot_matrix_x2

  !============================================================================

  function rot_matrix_y1(Angle) result(Rot_DD)

    real, intent(in) :: Angle
    real             :: Rot_DD(3,3)

    Rot_DD = rot_matrix_y(sin(Angle),cos(Angle))

  end function rot_matrix_y1

  !============================================================================

  function rot_matrix_y2(SinAngle,CosAngle) result(Rot_DD)

    real, intent(in) :: CosAngle, SinAngle
    real             :: Rot_DD(3,3)

    
    Rot_DD        =  0
    Rot_DD(z_,z_) =  CosAngle
    Rot_DD(z_,x_) = -SinAngle
    Rot_DD(x_,z_) =  SinAngle
    Rot_DD(x_,x_) =  CosAngle
    Rot_DD(y_,y_) =  1

  end function rot_matrix_y2

  !============================================================================

  function rot_matrix_z1(Angle) result(Rot_DD)

    real, intent(in) :: Angle
    real :: Rot_DD(3,3)

    rot_DD = rot_matrix_z(sin(Angle),cos(Angle))

  end function rot_matrix_z1

  !============================================================================

  function rot_matrix_z2(SinAngle,CosAngle) result(Rot_DD)

    real, intent(in) :: SinAngle,CosAngle
    real             :: Rot_DD(3,3)

    Rot_DD        =  0
    Rot_DD(x_,x_) =  CosAngle
    Rot_DD(x_,y_) = -SinAngle
    Rot_DD(y_,x_) =  SinAngle
    Rot_DD(y_,y_) =  CosAngle
    Rot_DD(z_,z_) =  1

  end function rot_matrix_z2

  !============================================================================

  function rot_xyz_sph1(Xyz_D) result(Rot_DD)

    real, intent(in) :: Xyz_D(3)
    real             :: Rot_DD(3,3)

    real :: SinTheta, CosTheta, SinPhi, CosPhi

    call xyz_to_dir(Xyz_D(1),Xyz_D(2),Xyz_D(3),SinTheta,CosTheta,SinPhi,CosPhi)

    Rot_DD = rot_xyz_sph(SinTheta,CosTheta,SinPhi,CosPhi)

  end function rot_xyz_sph1

  !============================================================================

  function rot_xyz_sph3(x,y,z) result(Rot_DD)

    real, intent(in) :: x,y,z
    real             :: Rot_DD(3,3)

    real :: SinTheta, CosTheta, SinPhi, CosPhi

    call xyz_to_dir(x,y,z,SinTheta,CosTheta,SinPhi,CosPhi)

    Rot_DD = rot_xyz_sph(SinTheta,CosTheta,SinPhi,CosPhi)

  end function rot_xyz_sph3

  !============================================================================

  function rot_xyz_sph2(Theta,Phi) result(Rot_DD)

    real, intent(in) :: Theta,Phi
    real             :: Rot_DD(3,3)

    Rot_DD = rot_xyz_sph4(sin(Theta),cos(Theta),sin(Phi),cos(Phi))

  end function rot_xyz_sph2

  !BOP -------------------------------------------------------------------
  !IROUTINE: rot_xyz_sph - Cartesian/spherical vector tranformation matrix
  !
  !INTERFACE:
  function rot_xyz_sph4(SinTheta,CosTheta,SinPhi,CosPhi) result(XyzSph_DD)

    !INPUT PARAMETERS:
    real, intent(in) :: SinTheta,CosTheta,SinPhi,CosPhi

    !RETURN VALUE:
    real             :: XyzSph_DD(3,3)

    !DESCRIPTION:
    !
    ! This is the fundamental routine that calculates the rotation
    ! matrix between Cartesian and spherical coordinates.
    ! The spherical direction of the location is given by the 4 
    ! trigonometric function arguments. 
    ! The matrix is obtained from 3 rotations:
    ! $$
    !      R = R_z(\theta-\pi/2) \cdot R_y(-\phi) \cdot R_x(-\pi/2)
    ! $$
    ! The resulting matrix is explicitly given here for sake of speed.
    !
    ! The subroutine should be called in one of the following ways:
    ! \begin{verbatim}
    !     XyzSph_DD = rot_xyz_sph(XyzPos_D)
    !
    !     XyzSph_DD = rot_xyz_sph(xPos,yPos,zPos)
    !
    !     XyzSph_DD = rot_xyz_sph(Theta,Phi)
    !
    !     XyzSph_DD = rot_xyz_sph(SinTheta,CosTheta,SinPhi,CosPhi)
    ! \end{verbatim}
    ! The 4 argument version requires the minimum amount of computations.
    ! The matrix can be used to convert back and forth like this:
    ! \begin{verbatim}
    !     Xyz_D = matmul(XyzSph_DD, Sph_D)
    !
    !     Sph_D = matmul(Xyz_D, XyzSph_DD)
    ! \end{verbatim}
    !EOP
    !----------------------------------------------------------------------
    !BOC
    XyzSph_DD = reshape( (/ &
         CosPhi*SinTheta, SinPhi*SinTheta,  CosTheta, &
         CosPhi*CosTheta, SinPhi*CosTheta, -SinTheta, &
         -SinPhi,         CosPhi,           0.0 /), &
         (/3,3/) )
    !EOC
  end function rot_xyz_sph4

  !BOP ========================================================================
  !IROUTINE: cross_product - return cross product of two vectors
  !INTERFACE:
  function cross_product11(a_D, b_D) result(c_D)
    !INPUT ARGUMENTS:
    real, intent(in) :: a_D(3), b_D(3)
    !RETURN VALUE:
    real             :: c_D(3)
    !EOP
    !-------------------------------------------------------------------------
    !BOC
    c_D(x_) = a_D(y_)*b_D(z_) - a_D(z_)*b_D(y_)
    c_D(y_) = a_D(z_)*b_D(x_) - a_D(x_)*b_D(z_)
    c_D(z_) = a_D(x_)*b_D(y_) - a_D(y_)*b_D(x_)
    !EOC
  end function cross_product11

  !============================================================================

  function cross_product13(a_D, bX, bY, bZ) result(c_D)
    real, intent(in) :: a_D(3), bX, bY, bZ
    real             :: c_D(3)
    !-------------------------------------------------------------------------
    c_D(x_) = a_D(y_)*bz - a_D(z_)*by
    c_D(y_) = a_D(z_)*bx - a_D(x_)*bz
    c_D(z_) = a_D(x_)*by - a_D(y_)*bx
  end function cross_product13

  !============================================================================

  function cross_product31(aX, aY, aZ, b_D) result(c_D)
    real, intent(in) :: aX, aY, aZ, b_D(3)
    real             :: c_D(3)
    !-------------------------------------------------------------------------
    c_D(x_) = ay*b_D(z_) - az*b_D(y_)
    c_D(y_) = az*b_D(x_) - ax*b_D(z_)
    c_D(z_) = ax*b_D(y_) - ay*b_D(x_)
  end function cross_product31

  !============================================================================

  function cross_product33(aX, aY, aZ, bX, bY, bZ) result(c_D)
    real, intent(in) :: aX, aY, aZ, bX, bY, bZ
    real             :: c_D(3)
    !-------------------------------------------------------------------------
    c_D(x_) = ay*bz - az*by
    c_D(y_) = az*bx - ax*bz
    c_D(z_) = ax*by - ay*bx
  end function cross_product33

  !============================================================================
  function inverse_matrix(a_DD, SingularLimit, DoIgnoreSingular) result(b_DD)

    ! Return the inverse of the 3x3 matrix a_DD
    ! The optional SingularLimit gives the smallest value for the determinant
    ! The optional DoIgnoreSingular determines what to do if the 
    ! determinant is less than SingularLimit.

    real, intent(in) :: a_DD(3,3)
    real, intent(in), optional :: SingularLimit
    logical, intent(in), optional :: DoIgnoreSingular

    real             :: b_DD(3,3)

    real    :: DetA, Limit
    logical :: DoIgnore

    character (len=*), parameter :: NameSub = NameMod//':inverse_matrix'
    !-------------------------------------------------------------------------
    Limit = 1.e-16
    if(present(SingularLimit)) Limit = SingularLimit
    DoIgnore = .false.
    if(present(DoIgnoreSingular)) DoIgnore = DoIgnoreSingular

    ! Invert the 3x3 matrix:
    b_DD(1,1)=a_DD(2,2)*a_DD(3,3)-a_DD(2,3)*a_DD(3,2)
    b_DD(2,1)=a_DD(2,3)*a_DD(3,1)-a_DD(2,1)*a_DD(3,3)
    b_DD(3,1)=a_DD(2,1)*a_DD(3,2)-a_DD(2,2)*a_DD(3,1)

    b_DD(1,2)=a_DD(1,3)*a_DD(3,2)-a_DD(1,2)*a_DD(3,3)
    b_DD(2,2)=a_DD(1,1)*a_DD(3,3)-a_DD(1,3)*a_DD(3,1)
    b_DD(3,2)=a_DD(1,2)*a_DD(3,1)-a_DD(1,1)*a_DD(3,2)

    b_DD(1,3)=a_DD(1,2)*a_DD(2,3)-a_DD(1,3)*a_DD(2,2)
    b_DD(2,3)=a_DD(1,3)*a_DD(2,1)-a_DD(1,1)*a_DD(2,3)
    b_DD(3,3)=a_DD(1,1)*a_DD(2,2)-a_DD(1,2)*a_DD(2,1)

    DetA= a_DD(1,1)*b_DD(1,1)+a_DD(1,2)*b_DD(2,1)+a_DD(1,3)*b_DD(3,1)

    if(abs(detA) > Limit)then
       b_DD = b_DD/DetA
    elseif(DoIgnore)then
       b_DD = cUnit_DD
    else
       write(*,*)'Error in ',NameSub,' for matrix:'
       call show_rot_matrix(a_DD)
       call CON_stop('Singular matrix in '//NameSub)
    end if

  end function inverse_matrix

  !============================================================================

  subroutine show_rot_matrix(Matrix_DD)

    real, intent(in) :: Matrix_DD(3,3)

    write(*,'(3(3f14.10,/))') transpose(Matrix_DD)

  end subroutine show_rot_matrix
  !BOP =======================================================================
  !IROUTINE: atan2_check - calculate atan2 with a check for 0.,0. arguments
  !INTERFACE:
  real function atan2_check(x,y)
    !INPUT ARGUMENTS:
    real, intent(in) :: x,y
    !EOP
    !BOC
    if(x==0 .and. y==0)then
       atan2_check = 0
    else
       atan2_check = atan2(x,y)
    end if
    if(atan2_check < 0.0)atan2_check = atan2_check + cTwoPi
    !EOC
  end function atan2_check

  !============================================================================

  subroutine test_coord_transform

    real, parameter      :: cTiny = 0.000001
    real, dimension(3)   :: Xyz_D,Sph_D,Xyz2_D
    real, dimension(3,3) :: XyzSph_DD

    Xyz_D = (/0.1, 0.2, 0.3/)
    write(*,'(a,3es16.8)')'Xyz_D=',Xyz_D

    call xyz_to_sph(Xyz_D,Sph_D)

    write(*,'(a,3es16.8)')'Sph_D=',Sph_D

    call sph_to_xyz(Sph_D,Xyz2_D)

    write(*,'(a,3es16.8)')'Xyz_D=',Xyz_D

    if(maxval(abs(Xyz_D-Xyz2_D)) > cTiny) &
         write(*,'(a)')'Error transforming xyz->sph->xyz'

    write(*,'(a,/,3(3f14.10,/))') 'rot_matrix_z(-Phi)='
    call show_rot_matrix(rot_matrix_z(-Sph_D(3)))

    Xyz_D = matmul(rot_matrix_z(-Sph_D(3)),Xyz_D)

    write(*,'(a,3es16.8)')'rot_matrix_z(-Phi).Xyz_D=',Xyz_D

    Xyz_D = matmul(Xyz_D,rot_matrix_y(Sph_D(2)))

    write(*,'(a,3es16.8)')'Xyz_D.rot_matrix_y(Theta)=',Xyz_D

    if(any(abs(Xyz_D(1:2)) > cTiny)) &
         write(*,'(a)')'Error rotating Xyz_D into Z axis'

    if(abs(Xyz_D(3)-Sph_D(1)) > cTiny) &
         write(*,'(a)')'Error rotating Xyz_D, length changed'

    Xyz_D = (/0.001, -0.4, 0.35353/)
    write(*,'(a,3es16.8)')'Original Xyz=',Xyz_D
    Xyz2_D = matmul(rot_matrix_x(1.),Xyz_D)
    Xyz2_D = matmul(rot_matrix_y(2.),Xyz2_D)
    Xyz2_D = matmul(rot_matrix_z(3.),Xyz2_D)
    Xyz2_D = matmul(rot_matrix_z(-3.),Xyz2_D)
    Xyz2_D = matmul(rot_matrix_y(-2.),Xyz2_D)
    Xyz2_D = matmul(rot_matrix_x(-1.),Xyz2_D)
    write(*,'(a,3es16.8)')'Rotated  Xyz=',Xyz2_D

    if(maxval(abs(Xyz_D-Xyz2_D)) > cTiny) &
         write(*,'(a)')'Error rotating back and forth'

    write(*,*)
!    Xyz_D = (/1.0, 0.0, 0.0/)
!    Xyz_D = (/8.0, 6.0, 0.0/)
!    Xyz_D = (/8.0, 0.0, 6.0/)
    Xyz_D = (/8.0, 0.1, 6.0/)
    write(*,'(a,3es16.8)')'Cartesian position=',Xyz_D
    XyzSph_DD = rot_xyz_sph(Xyz_D)
    write(*,'(a)')'XyzSph_DD'; call show_rot_matrix(XyzSph_DD)
    Sph_D = (/1.0, 0.0, 0.0/)
    write(*,'(a,3es16.8)')'Spherical vector  =',Sph_D
    write(*,'(a,3es16.8)')'Cartesian vector  =',matmul(XyzSph_DD,Sph_D)
    Sph_D = (/0.0, 1.0, 0.0/)
    write(*,'(a,3es16.8)')'Spherical vector  =',Sph_D
    write(*,'(a,3es16.8)')'Cartesian vector  =',matmul(XyzSph_DD,Sph_D)
    Sph_D = (/0.0, 0.0, 1.0/)
    write(*,'(a,3es16.8)')'Spherical vector  =',Sph_D
    write(*,'(a,3es16.8)')'Cartesian vector  =',matmul(XyzSph_DD,Sph_D)

  end subroutine test_coord_transform

end module ModCoordTransform
