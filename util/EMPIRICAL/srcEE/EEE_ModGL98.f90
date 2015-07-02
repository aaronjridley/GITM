!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==============================================================================
module EEE_ModGL98
  use EEE_ModCommonVariables
  implicit none
  save

  private

  public :: set_parameters_GL98
  public :: get_GL98_fluxrope
  public :: adjust_GL98_fluxrope

  real :: cme_a=0.0, cme_r1=0.0, cme_r0=0.0, cme_a1=0.0, cme_alpha=0.0
  real :: cme_rho1=0.0, cme_rho2=0.0

  real :: ModulationRho=0.0, ModulationP=0.0

contains

  !============================================================================

  subroutine set_parameters_GL98(NameCommand)
    use ModReadParam, ONLY: read_var
    character(len=*), intent(in):: NameCommand

    character(len=*), parameter:: NameSub = 'set_parameters_GL98'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case("#GL98FLUXROPE")
       call read_var('UseFluxRope',     DoAddFluxRope)
       call read_var('cme_a',           cme_a)
       call read_var('cme_r1',          cme_r1)
       call read_var('cme_r0',          cme_r0)
       call read_var('cme_a1',          cme_a1)
       call read_var('cme_alpha',       cme_alpha)
       call read_var('cme_rho1',        cme_rho1)
       call read_var('cme_rho2',        cme_rho2)
       call read_var('ModulationRho',   ModulationRho)
       call read_var('ModulationP',     ModulationP)
       call read_var('LongitudeCme',   LongitudeCme)
       call read_var('LatitudeCme',    LatitudeCme)
       call read_var('OrientationCme', OrientationCme)
    case("#CME")
       call read_var('Stretch',     cme_a)
       call read_var('Distance',    cme_r1)
       call read_var('Radius',      cme_r0)
       call read_var('Bstrength',   cme_a1)
       call read_var('Density',     cme_rho1)
       call read_var('ModulationRho',   ModulationRho)
       call read_var('ModulationP',     ModulationP)
    case default
       call CON_stop(NameSub//' unknown NameCommand='//NameCommand)
    end select

  end subroutine set_parameters_GL98

  !============================================================================

  subroutine get_GL98_fluxrope(R_GL98_D,rho_GL98,p_GL98,B_GL98_D)
    !--------------------------------------------------------------------------
    ! Definition of Parameters used for the initial state
    !   cme_a    = contraction distance as in   r --> r -a
    !   cme_r1   = distance of flux rope from sun center = 1.2
    !   cme_r0   = radius of flux rope
    !   cme_a1   = constant for setting pressure in flux rope
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////
    !=====================================================================
    !\
    ! Calculates magnetic field, pressure and density for a coronal flux
    ! rope capable of self-similar expansion and erupting in a CME.
    ! The analytical solution is taken from Gibson and Low 
    ! Astrophysical Journal, Vol 493, p. 460.
    !
    ! Written by Chip Manchester Jan 18 2001
    ! Rewritten by Chip Nov 29 2001 for flux rope injection
    !/
    !   Bug fixes
    !   March 18       dpres_1dr1 = cme_a1*(dA2dr*dr2dr1 + dA2dth*dth2dr1)
    !   change to..... dpres_1dr1 = a1scl*(dA2dr*dr2dr1 + dA2dth*dth2dr1)
    !   without above fix, same as used for runs 12, 13, 14
    !   Feb  07, 2002 Br1 is changed to Br1**2 in density calc thanks to Ilia
    !   Feb  12, 2002 expression for ga0r0 is fixed     
    !   Feb  22, 2002 derivative of B2*2/dr1 is fixed
    !   Feb  22, 2002 angles in 2nd coordinate system are redefined
    !
    !------------------------------------------------------------------------
    use ModNumConst,       ONLY: cPi, cDegToRad
    use ModCoordTransform, ONLY: rot_matrix_y, rot_matrix_z
    implicit none
    real, dimension(3), intent(in) :: R_GL98_D
    real, dimension(3), intent(out) :: B_GL98_D
    real, intent(out) :: rho_GL98,p_GL98
    !\
    ! User declared local variables go here::
    !/
    real :: x,y,z,  x_1,y_1,z_1,  x_2,y_2,z_2
    real :: r,cos_theta,sin_theta,cos_phi,sin_phi
    real :: r_1,cos_theta1,sin_theta1,cos_phi1,sin_phi1,lambda
    real :: r_2,cos_theta2,sin_theta2,cos_phi2,sin_phi2
    real :: dr2dr1,dth2dr1
    real :: Br,  Btheta, Bphi
    real :: Br1, Btheta1, Bphi1 
    real :: Br2, Btheta2, Bphi2
    real :: Bx_1, By_1, Bz_1
    real :: Bx_2, By_2, Bz_2
    real :: dBr2dr1, dBtheta2dr1,dBphi2dr1
    real :: dBr2dr2, dBth2dr2, dBphi2dr2
    real :: dBr2dth2, dBth2dth2, dBphi2dth2
    real :: A2,dA2dr, dA2dth, d2A2dr2, d2A2drdth, d2A2dth2
    real :: pres_1,dpres_1dr1,F_grav,alpha0,ga0r0,delta 
    real, dimension(3) :: R1_GL98_D,B1_GL98_D
    real, dimension(3,3), save :: RotateGL98_DD
    logical, save :: DoFirst_GL=.true.

    real :: a1scl,rho1scl,rho2scl,SSPscl
    !------------------------------------------------------------------------
    if (DoFirst_GL) then
       DoFirst_GL=.false.
       !\
       ! Construct the rotational matrix RotateGL98_DD,
       !/
       RotateGL98_DD  = matmul( &
            rot_matrix_z(OrientationCme*cDegToRad),&
            rot_matrix_y((LatitudeCme-90)*cDegToRad))
       RotateGL98_DD = matmul(RotateGL98_DD, &
            rot_matrix_z(-LongitudeCme*cDegToRad))

       if(iProc==0)then
          write(*,*) prefix
          write(*,*) prefix, &
               '>>>>>>>>>>>>>>>>>>>                  <<<<<<<<<<<<<<<<<<<<<'
          write(*,*) prefix
          write(*,*) prefix, &
               '            Gibson and Low CME is initiated!!!'
          write(*,*) prefix
          write(*,*) prefix, &
               '>>>>>>>>>>>>>>>>>>>                  <<<<<<<<<<<<<<<<<<<<<'
          write(*,*) prefix
          write(*,*) prefix, 'cme_a  = ',cme_a, '[rSun]'
          write(*,*) prefix, 'cme_r1 = ',cme_r1,'[rSun]'
          write(*,*) prefix, 'cme_r0 = ',cme_r0,'[rSun]'
          write(*,*) prefix, 'cme_a1 = ',cme_a1,'[Gauss]'
          write(*,*) prefix, 'cme_rho1 = ',cme_rho1,'[kg/m^3]'
          write(*,*) prefix, 'cme_rho2 = ',cme_rho2,'[kg/m^3]'
          write(*,*) prefix, 'ModulationRho = ',ModulationRho
          write(*,*) prefix, 'ModulationP   = ',ModulationP
          write(*,*) prefix, 'LongitudeCme = ',LongitudeCme,'[degrees]'
          write(*,*) prefix, 'LatitudeCme = ',LatitudeCme,'[degrees]'
          write(*,*) prefix, 'OrientationCme = ',OrientationCme,'[degrees]'
          write(*,*) prefix
       end if
    end if

    a1scl   = cme_a1*Io2No_V(UnitB_)
    rho1scl = cme_rho1*Si2No_V(UnitRho_)
    rho2scl = cme_rho2*Si2No_V(UnitRho_)
    SSPscl  = 1.0

    delta = 0.1
    !\
    ! Compute R1_GL98_D::
    !/
    R1_GL98_D = matmul(RotateGL98_DD, R_GL98_D)
    !\
    ! CALCULATE CELL CENTER FOR GLOBAL CARTESIAN COORDINATES 
    !/
    x = R1_GL98_D(x_)
    y = R1_GL98_D(y_)
    z = R1_GL98_D(z_)
    !\
    ! CALCULATE CELL CENTER FOR GLOBAL SPHERICAL COORDINATES 
    !/
    r = sqrt(x**2 + y**2 + z**2)
    cos_theta = z/r
    sin_theta = sqrt(x**2 + y**2)/r
    cos_phi   = x/sqrt(x**2 + y**2)
    sin_phi   = y/sqrt(x**2 + y**2)
    if (r <= delta) then 
       r = delta
       x = delta*sin_theta*cos_phi
       y = delta*sin_theta*sin_phi
       z = delta*cos_theta
    end if
    !\
    ! CALCULATE CELL CENTER FOR TRANSFORMED SPHERICAL COORDINATES 
    ! stretching transformation of variables r --> r - a
    !/
    lambda = r + cme_a
    r_1    = lambda
    cos_theta1 = cos_theta
    sin_theta1 = sin_theta 
    cos_phi1   = cos_phi
    sin_phi1   = sin_phi
    !\
    ! CALCULATE CELL CENTER FOR TRANSFORMED CARTESIAN COORDINATES 
    !/
    x_1 = lambda*sin_theta1*cos_phi1
    y_1 = lambda*sin_theta1*sin_phi1
    z_1 = lambda*cos_theta1
    !---------------------------FLUX ROPE REGION---------------------
    ! CALCULATE CELL CENTER CARTESIAN COORDINATES for CME FLUX ROPE
    ! stretching transformation r = r --> r - a
    !----------------------------------------------------------------
    x_2 = x_1
    y_2 = y_1 
    z_2 = z_1 - cme_r1 
    !\
    ! CALCULATE CELL CENTER SPHERICAL COORDINATES for CME FLUX ROPE
    !/
    r_2 = sqrt(x_2**2 + y_2**2 + z_2**2)
    cos_theta2 = x_2/r_2
    sin_theta2 = sqrt(z_2**2 + y_2**2)/r_2
    cos_phi2   = y_2/sqrt(z_2**2 + y_2**2)
    sin_phi2   = z_2/sqrt(z_2**2 + y_2**2)
    if (r_2 <= delta) then
       r_2 = delta
       y_2 = delta*sin_theta2*cos_phi2
       z_2 = delta*sin_theta2*sin_phi2
       x_2 = delta*cos_theta2
    end if
    alpha0 = 5.763854/cme_r0
    ga0r0 = sin(alpha0*cme_r0)/(alpha0*cme_r0) - cos(alpha0*cme_r0)
    A2 = (4.0*cPi*a1scl/alpha0**2)*((cme_r0**2/ga0r0) &
         *(sin(alpha0*r_2)/(alpha0*r_2) - cos(alpha0*r_2)) - r_2**2) &
         *sin_theta2**2
    dA2dr = ((4.0*cPi*a1scl/alpha0**2)*((cme_r0**2/ga0r0) &
         *(cos(alpha0*r_2)/r_2 - sin(alpha0*r_2)/(alpha0*r_2**2) &
         + alpha0*sin(alpha0*r_2)) - 2.0*r_2))*sin_theta2**2
    dA2dth = (8.0*cPi*a1scl/alpha0**2)*((cme_r0**2/ga0r0) &
         *(sin(alpha0*r_2)/(alpha0*r_2) - cos(alpha0*r_2)) - r_2**2) &
         *sin_theta2*cos_theta2 
    d2A2dr2 = (4.0*cPi*a1scl/alpha0**2)*sin_theta2**2 &
         *( (cme_r0**2/ga0r0)*(2.0*sin(alpha0*r_2)/(alpha0*r_2**3) &
         - 2.0*cos(alpha0*r_2)/(r_2**2) - alpha0*sin(alpha0*r_2)/r_2 &
         + (alpha0**2)*cos(alpha0*r_2)) - 2.0)  
    d2A2drdth = (8.0*cPi*a1scl/alpha0**2)*sin_theta2*cos_theta2 &
         *((cme_r0**2/ga0r0)*(cos(alpha0*r_2)/r_2 &
         - sin(alpha0*r_2)/(alpha0*r_2**2) &
         + alpha0*sin(alpha0*r_2)) - 2.0*r_2) 
    d2A2dth2 = (8.0*cPi*a1scl/alpha0**2)*((cme_r0**2/ga0r0) &
         *(sin(alpha0*r_2)/(alpha0*r_2) - cos(alpha0*r_2)) - r_2**2) &
         *(cos_theta2**2 - sin_theta2**2)
    dr2dr1  =  (x_2/r_2)*sin_theta1*cos_phi1 &
         + (y_2/r_2)*sin_theta1*sin_phi1 &
         + (z_2/r_2)*cos_theta1
    dth2dr1 =  (1.0/r_2)*(-sin_theta2*sin_theta1*cos_phi1 &
         + cos_theta2*cos_phi2*sin_theta1*sin_phi1 &
         + cos_theta2*sin_phi2*cos_theta1)
    !\
    ! Derivatives of field components in flux rope spherical coordinates
    !/
    dBr2dr2 = -2.0*dA2dth/(sin_theta2*r_2**3) &
         + d2A2drdth/(sin_theta2*r_2**2)
    dBr2dth2 = -cos_theta2*dA2dth/(r_2**2*sin_theta2**2) & 
         + d2A2dth2/(sin_theta2*r_2**2)
    dBth2dr2 = dA2dr/(sin_theta2*r_2**2) &
         - d2A2dr2/(r_2 *sin_theta2)
    dBth2dth2 = cos_theta2*dA2dr/(r_2*sin_theta2**2) &
         - d2A2drdth/(r_2*sin_theta2)
    dBphi2dr2 = alpha0*dA2dr/(r_2*sin_theta2) &
         - alpha0*A2 /(sin_theta2*r_2**2)
    dBphi2dth2 = alpha0*dA2dth/(r_2*sin_theta2) &
         - alpha0*cos_theta2*A2/(r_2*sin_theta2**2)
    !\
    ! Total derivative of the flux rope field components in terms of `r1'
    !/
    dBr2dr1     = dBr2dr2  *dr2dr1 + dBr2dth2  *dth2dr1
    dBtheta2dr1 = dBth2dr2 *dr2dr1 + dBth2dth2 *dth2dr1
    dBphi2dr1   = dBphi2dr2*dr2dr1 + dBphi2dth2*dth2dr1
    !\
    ! Magnetic field components in the flux rope spherical coordinates
    !/
    Br2     = dA2dth/(sin_theta2*r_2**2)
    Btheta2 = -dA2dr/(sin_theta2*r_2)
    Bphi2   = alpha0*A2/(sin_theta2*r_2)
    !\
    ! Magnetic field components in the second cartesian coordinates
    ! X-COMPONENT OF MAGNETIC FIELD
    !/
    Bx_2 = Br2*cos_theta2 &
         - Btheta2*sin_theta2
    ! Y-COMPONENT OF MAGNETIC FIELD
    By_2 = Br2*sin_theta2*cos_phi2 &
         + Btheta2*cos_theta2*cos_phi2 &
         - Bphi2*sin_phi2
    ! Z-COMPONENT OF MAGNETIC FIELD
    Bz_2 = Br2*sin_theta2*sin_phi2 &
         + Btheta2*cos_theta2*sin_phi2 &
         + Bphi2*cos_phi2
    !\
    ! Define the magnetic field in the global cartesian coordinates
    ! INSIDE THE MAGNETIC FLUX ROPE REGION
    !/
    if (sqrt(x_2**2 + y_2**2 + z_2**2) <= cme_r0) then
       Bx_1 = Bx_2 
       By_1 = By_2 
       Bz_1 = Bz_2 
       !\
       ! Magnetic field components in global sperical coordinates
       !/
       Br1 = Bx_1*sin_theta1*cos_phi1 &
            + By_1*sin_theta1*sin_phi1 &
            + Bz_1*cos_theta1
       Btheta1 = Bx_1*cos_theta1*cos_phi1 &
            + By_1*cos_theta1*sin_phi1 &
            - Bz_1*sin_theta1
       Bphi1 = -Bx_1*sin_phi1 &
            + By_1*cos_phi1
       !\
       ! Compute kinetic gas pressure
       !/
       pres_1     = inv_g*rho1scl*SSPscl**2 + a1scl*A2
       dpres_1dr1 = a1scl*(dA2dr*dr2dr1 + dA2dth*dth2dr1)
       !\
       ! MAGNETIC FIELD transformed with stretching transformation
       !/
       Br     = Br1    *(lambda/r)**2
       Btheta = Btheta1*(lambda/r)
       Bphi   = Bphi1  *(lambda/r)
       !\
       ! Magnetic field components in global cartesian coordinates
       !/
       B1_GL98_D(x_) = Br*sin_theta1*cos_phi1 &
            + Btheta*cos_theta1*cos_phi1 &
            - Bphi*sin_phi1
       !\
       ! Y-COMPONENT OF MAGNETIC FIELD:: (x,y,z)_BATSRUS -> (x,y,z)
       !/
       B1_GL98_D(y_) = Br*sin_theta1*sin_phi1 &
            + Btheta*cos_theta1*sin_phi1 &
            + Bphi*cos_phi1
       !\
       ! Z-COMPONENT OF MAGNETIC FIELD:: (x,y,z)_BATSRUS -> (x,y,z)
       !/
       B1_GL98_D(z_) = Br*cos_theta1 - Btheta*sin_theta1
       !\
       ! Transform back to the original coordinates
       ! given by R_GL98_D:: 
       !/
       B_GL98_D  = matmul(B1_GL98_D, RotateGL98_DD)
       !\
       ! PLASMA DENSITY with stretching transformation
       !/
       F_grav   = (abs(Gbody)/(r**2) + cme_alpha*r)
       rho_GL98 = (1.0/ F_grav)*(((lambda/r)**2) &
            *((lambda/r)**2 - 1.0)*(dpres_1dr1 + (1.0/(4.0*cPi)) &
            *(Br2*dBr2dr1 + Btheta2*dBtheta2dr1 + Bphi2*dBphi2dr1)) &
            + 2.0*lambda*cme_a*pres_1/r**3 + cme_a*lambda &
            /(4.0*cPi*r**3)*(1.0 - 2.0*(lambda/r)**2)*Br1**2 &
            + ((lambda/r)**2)*((cme_a/r)**2 +2.0*cme_a/r)*(Btheta1**2 &
            + Bphi1**2)/(4.0*cPi*lambda))
       !\
       ! Add background density
       !/
       rho_GL98 = rho_GL98 + rho2scl/r**3
       !\
       ! PLASMA PRESSURE with contraction transformation
       !/
       p_GL98 = pres_1*(lambda/r)**2 &
            - (1.0/(8.0*cPi))*((lambda/r)**2)*((lambda/r)**2 - 1.0)*Br1**2 
       !\
       ! Add background pressure
       !/
       p_GL98 = p_GL98 + abs(Gbody)*rho2scl/(4.0*r**4)

       rho_GL98 = rho_GL98*No2Si_V(UnitRho_)
       p_GL98 = p_GL98*No2Si_V(UnitP_)
       B_GL98_D = B_GL98_D*No2Si_V(UnitB_)

    else
       B_GL98_D = 0.0; rho_GL98 = 0.0; p_GL98 = 0.0
    endif

  end subroutine get_GL98_fluxrope

  !============================================================================

  subroutine adjust_GL98_fluxrope(Rho,p)

    real, intent(inout) :: Rho,p
    !--------------------------------------------------------------------------

    !\
    ! Add just `ModulationRho' times of the CME mass
    ! to the mass density::
    !/
    if(ModulationRho*Rho <= 0.0)then
       Rho = Rho
    else
       Rho = ModulationRho*Rho
    end if

    !\
    ! Add just `ModulationP' times of the CME pressure
    ! to the kinetic pressure::
    !/
    if(ModulationP*p <= 0.0) then
       p = p
    else
       p = ModulationP*p
    end if

  end subroutine adjust_GL98_fluxrope

end module EEE_ModGL98
