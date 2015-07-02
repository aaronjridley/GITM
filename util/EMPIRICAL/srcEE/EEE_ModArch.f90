!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==============================================================================
module EEE_ModArch
  use EEE_ModCommonVariables
  implicit none
  save

  private

  public :: set_parameters_arch, get_arch

  real    :: DipoleStrengthSi,rDipole,AngleBetweenDipoles
  integer :: nDipole

contains

  !============================================================================

  subroutine set_parameters_arch(NameCommand)

    use ModReadParam, ONLY: read_var

    character(len=*), intent(in):: NameCommand

    character(len=*), parameter:: NameSub = 'set_parameters_arch'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case("#ARCH")
       call read_var('UseArch',          UseArch)
       call read_var('DipoleStrengthSi', DipoleStrengthSi)
       call read_var('rDipole',          rDipole)
       call read_var('nDipole', nDipole)
       if(nDipole > 1) &
            call read_var('AngleBetweenDipoles', AngleBetweenDipoles)
       call read_var('LongitudeCme',   LongitudeCme)
       call read_var('LatitudeCme',    LatitudeCme)
       call read_var('OrientationCme', OrientationCme)
    case("#CME")
       call read_var('DipoleStrengthSi', DipoleStrengthSi)
       call read_var('rDipole',          rDipole)
       call read_var('nDipole',          nDipole)
       if(nDipole > 1) &
            call read_var('AngleBetweenDipoles', AngleBetweenDipoles)
    case default
       call CON_stop(NameSub//' unknown NameCommand='//NameCommand)
    end select

    if(rDipole>=1.0 .or. rDipole<0.0)then
       call CON_stop('rDipole in #ARCH should be inside the Sun')
    end if
    if(nDipole < 1)then
       call CON_stop('nDipole need to be larger than 0 in #ARCH')
    end if

  end subroutine set_parameters_arch

  !============================================================================

  subroutine get_arch(x0_D,B0_D)
    use EEE_ModCommonVariables
    use ModCoordTransform, ONLY: rot_matrix_x,rot_matrix_y,rot_matrix_z
    implicit none

    real, intent(in) :: x0_D(3)
    real, intent(out) :: B0_D(3)

    logical, save :: DoFirst=.true.
    real, allocatable, dimension(:,:), save :: xDip_DI
    real, save :: Bdp
    real, dimension(3,3), save :: Rotate_DD
    real :: Phi,PhiFirstDipole
    real :: r,r_inv,r2_inv,r3_inv,x_D(3),xShift_D(3)
    real :: Dp,Bdp_D(3),B_D(3)
    integer :: iDipole
    !--------------------------------------------------------------------------
    if(DoFirst)then
       DoFirst = .false.

       allocate(xDip_DI(3,nDipole))

       PhiFirstDipole = -0.5*AngleBetweenDipoles*(nDipole - 1)
       do iDipole=1,nDipole
          Phi = PhiFirstDipole + AngleBetweenDipoles*(iDipole - 1)

          xDip_DI(x_,iDipole) = rDipole*cos(Phi*cDegToRad)
          xDip_DI(y_,iDipole) = rDipole*sin(Phi*cDegToRad)
          xDip_DI(z_,iDipole) = 0.0

          Bdp  = DipoleStrengthSi*Si2No_V(UnitB_)
       end do

       ! Orientation rotates clockwise
       Rotate_DD = matmul(rot_matrix_x(OrientationCme*cDegToRad), &
            rot_matrix_y(LatitudeCme*cDegToRad))
       Rotate_DD = matmul(Rotate_DD,rot_matrix_z(-LongitudeCme*cDegToRad))
    end if


    x_D = matmul(Rotate_DD,x0_D)

    B_D = 0.0
    do iDipole=1,nDipole
       xShift_D = x_D - xDip_DI(:,iDipole)

       r = max(sqrt(dot_product(xShift_D,xShift_D)),cTolerance)
       r_inv = 1.0/r
       r2_inv = r_inv*r_inv
       r3_inv = r_inv*r2_inv

       Bdp_D = Bdp*(/ 0.0, 0.0, 1.0 /)

       Dp = dot_product(Bdp_D,xShift_D)*3.0*r2_inv

       B_D = B_D + (Dp*xShift_D-Bdp_D)*r3_inv
    end do


    B0_D = matmul(B_D,Rotate_DD)

    B0_D = B0_D*No2Si_V(UnitB_)

  end subroutine get_arch

end module EEE_ModArch
