!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==============================================================================
module EEE_ModShearFlow
  implicit none
  save

  private

  public :: set_parameters_shearflow
  public :: get_shearflow

  logical, public :: UseShearFlow=.false.
  real :: FlowAmplitude, FlowWidthAngle, MaxBrActiveRegion
  real :: StartTime, StopTime, RampUpTime, RampDownTime

  real :: xFlow_D(3), MaxBr

contains

  !============================================================================

  subroutine set_parameters_shearflow(NameCommand)

    use ModReadParam, ONLY: read_var
    use EEE_ModCommonVariables, ONLY: LongitudeCme, LatitudeCme

    character(len=*), intent(in):: NameCommand

    character(len=*), parameter:: NameSub = 'set_parameters_shearflow'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case("#SHEARFLOW")
       call read_var('UseShearFlow',      UseShearFlow)
       call read_var('FlowAmplitude',     FlowAmplitude)
       call read_var('FlowWidthAngle',    FlowWidthAngle)
       call read_var('MaxBrActiveRegion', MaxBrActiveRegion)
       call read_var('StartTime',         StartTime)
       call read_var('StopTime',          StopTime)
       call read_var('RampUpTime',        RampUpTime)
       call read_var('RampDownTime',      RampDownTime)
       call read_var('LongitudeCme',      LongitudeCme)
       call read_var('LatitudeCme',       LatitudeCme)
    case("#CME")
       call read_var('FlowAmplitude',     FlowAmplitude)
       call read_var('FlowWidthAngle',    FlowWidthAngle)
       call read_var('MaxBrActiveRegion', MaxBrActiveRegion)
       call read_var('StartTime',         StartTime)
       call read_var('StopTime',          StopTime)
       call read_var('RampUpTime',        RampUpTime)
       call read_var('RampDownTime',      RampDownTime)
    case default
       call CON_stop(NameSub//' unknown NameCommand='//NameCommand)
    end select

    if(RampUpTime<=0.0) call CON_stop(NameSub// &
         ': RampUpTime in #SHEARFLOW should be greater than zero')

    if(RampDownTime<=0.0) call CON_stop(NameSub// &
         ': RampDownTime in #SHEARFLOW should be greater than zero')

  end subroutine set_parameters_shearflow

  !============================================================================

  subroutine get_shearflow(x_D,Time,U_D,iteration_number)

    !DESCRIPTION:
    ! Boundary shear flow that conserves the radial magnetic flux

    use EEE_ModCommonVariables
    use EEE_ModGetB0,   ONLY: EEE_get_B0
    use ModMagnetogram, ONLY: get_magnetogram_field
    use ModNumConst,    ONLY: cTolerance, cDegToRad
    implicit none

    real, intent(in)  :: x_D(3), Time
    real, intent(out) :: U_D(3)
    integer, intent(in) :: iteration_number

    real :: TimeProfile
    real :: R, Theta, Phi, Xy, SinTheta, CosTheta, SinPhi, CosPhi, UnitR_D(3)
    real :: dTheta, dPhi
    real :: B_D(3), B0_D(3), FullBr, FullBrL, FullBrR
    real :: ShearProfileL, ShearProfileR
    real :: UTheta, UPhi

    logical, save :: DoFirst=.true.
    !--------------------------------------------------------------------------
    if(DoFirst)then
       DoFirst = .false.
       xFlow_D(1) = cos(LatitudeCme*cDegToRad)*cos(LongitudeCme*cDegToRad)
       xFlow_D(2) = cos(LatitudeCme*cDegToRad)*sin(LongitudeCme*cDegToRad)
       xFlow_D(3) = sin(LatitudeCme*cDegToRad)

       FlowWidthAngle = FlowWidthAngle*cDegToRad
       MaxBr = abs(MaxBrActiveRegion)*Si2No_V(UnitB_)
    end if

    if(Time < StartTime .or. Time > StopTime)then
       U_D = 0.0
       return
    else
       if(Time > StopTime - RampDownTime)then
          TimeProfile = (StopTime - Time)/RampDownTime
       else
          TimeProfile = min((Time - StartTime)/RampUpTime,1.0)
       end if
    end if
    !!! TimeProfile = 1.0

    R        = sqrt(sum(x_D**2))
    Theta    = acos(x_D(3)/R)
    Phi      = atan2(x_D(2),x_D(1))
    Xy       = max(sqrt(x_D(1)**2+x_D(2)**2),cTolerance)
    SinTheta = Xy/R
    CosTheta = x_D(3)/R
    SinPhi   = x_D(2)/Xy
    CosPhi   = x_D(1)/Xy
    UnitR_D  = x_D/R

    call get_magnetogram_field(UnitR_D(1),UnitR_D(2),UnitR_D(3),B0_D)
    call EEE_get_B0(UnitR_D,B_D)
    B0_D = (B0_D + B_D)*Si2No_V(UnitB_)

    FullBr = dot_product(UnitR_D,B0_D)

    dTheta = cTiny
    dPhi   = cTiny

    if(abs(FullBr)<cTolerance)then
       ! No flow at polarity inversion lines
       U_D = 0.0
       return
    else
       ShearProfileL = shear_profile(R,Theta,Phi-0.5*dPhi,Time,FullBrL)
       ShearProfileR = shear_profile(R,Theta,Phi+0.5*dPhi,Time,FullBrR)
       if(FullBrL*FullBrR<0.0 .and. abs(FullBrL)>cTolerance &
            .and. abs(FullBrR)>cTolerance)then
          UTheta = 0.0
       else
          UTheta = 1.0/FullBr &
               *(ShearProfileR-ShearProfileL)/(R*SinTheta*dPhi)
       end if

       ShearProfileL = shear_profile(R,Theta-0.5*dTheta,Phi,Time,FullBrL)
       ShearProfileR = shear_profile(R,Theta+0.5*dTheta,Phi,Time,FullBrR)
       if(FullBrL*FullBrR<0.0 .and. abs(FullBrL)>cTolerance &
            .and. abs(FullBrR)>cTolerance)then
          UPhi = 0.0
       else
          UPhi = -1.0/FullBr &
               *(ShearProfileR-ShearProfileL)/(R*dTheta)
       end if
    end if

    U_D(1) = UTheta*CosTheta*CosPhi - UPhi*SinPhi
    U_D(2) = UTheta*CosTheta*SinPhi + UPhi*CosPhi
    U_D(3) =-UTheta*SinTheta

    if(iteration_number == -1)then
       U_D = FlowAmplitude*U_D
    else
       U_D = FlowAmplitude*U_D*TimeProfile
    end if

  end subroutine get_shearflow

  !============================================================================

  real function shear_profile(R,Theta,Phi,Time,FullBr)
    use EEE_ModCommonVariables, ONLY: Si2No_V, UnitB_
    use EEE_ModGetB0,   ONLY: EEE_get_B0
    use ModMagnetogram, ONLY: get_magnetogram_field
    implicit none

    real, intent(in)  :: R, Theta, Phi, Time
    real, intent(out) :: FullBr

    real :: Xy, x_D(3), UnitR_D(3)
    real :: B0_D(3), B_D(3), BrRatio
    real :: Del_D(3), Angle, Mask
    !--------------------------------------------------------------------------
    Xy = R*sin(Theta)
    x_D(1) = Xy*cos(Phi)
    x_D(2) = Xy*sin(Phi)
    x_D(3) = R*cos(Theta)
    UnitR_D = x_D/R

    call get_magnetogram_field(UnitR_D(1),UnitR_D(2),UnitR_D(3),B0_D)
    call EEE_get_B0(UnitR_D,B_D)
    B0_D = (B0_D + B_D)*Si2No_V(UnitB_)

    FullBr = dot_product(UnitR_D,B0_D)
    BrRatio = FullBr/MaxBr

    Del_D = xFlow_D - UnitR_D
    Angle = 2.0*asin(0.5*sqrt(dot_product(Del_D,Del_D)))
    Mask = exp(-(Angle/FlowWidthAngle)**2)

    shear_profile = FullBr**3*exp(-BrRatio**2)*Mask

  end function shear_profile

end module EEE_ModShearFlow
