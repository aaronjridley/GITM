!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==============================================================================
Module EEE_ModMain

  use EEE_ModGetB0, ONLY: EEE_get_B0

  implicit none
  save

  private

  public :: EEE_initialize
  public :: EEE_set_parameters
  public :: EEE_get_state_init
  public :: EEE_get_state_BC
  public :: EEE_get_B0

contains
  !============================================================================
  subroutine EEE_initialize(BodyNDim,BodyTDim,gamma)
    use EEE_ModCommonVariables
    implicit none

    real, intent(in) :: BodyNDim,BodyTDim,gamma

    integer :: iComm, iError
    !--------------------------------------------------------------------------

    iComm = MPI_COMM_WORLD
    call MPI_COMM_RANK(iComm,iProc,iError)

    g = gamma
    inv_g = 1.0/g
    gm1 = g - 1.0
    inv_gm1 = 1.0/(g - 1.0)

    ! assume MassIon_I(1) = 1.0
    No2Si_V(UnitX_)   = rSun
    No2Si_V(UnitU_)   = sqrt(g*cBoltzmann*BodyTDim/cProtonMass)
    No2Si_V(UnitRho_) = 1000000*cProtonMass*BodyNDim

    !\
    ! Set other normalizing SI variables from the independent ones.
    !
    ! For sake of convenience
    !  units of B are chosen to satisfy v_A = B/sqrt(rho)       (mu = 1)
    !  units of n are chosen to satisfy  n  = rho/(ionmass/amu) (mp = 1)
    !  units of T are chosen to satisfy  T  = p/n               (kBoltzmann = 1)
    !
    ! Note that No2Si_V(UnitN_) is NOT EQUAL TO 1/No2Si_V(UnitX_)^3 !!!
    !/
    No2Si_V(UnitT_)          = No2Si_V(UnitX_)/No2Si_V(UnitU_)         ! s
    No2Si_V(UnitN_)          = No2Si_V(UnitRho_)/cProtonMass           ! #/m^3
    No2Si_V(UnitP_)          = No2Si_V(UnitRho_)*No2Si_V(UnitU_)**2    ! Pa
    No2Si_V(UnitB_)          = No2Si_V(UnitU_) &
         *sqrt(cMu*No2Si_V(UnitRho_))                                  ! T
    No2Si_V(UnitRhoU_)       = No2Si_V(UnitRho_)*No2Si_V(UnitU_)       ! kg/m^2/s
    No2Si_V(UnitEnergyDens_) = No2Si_V(UnitP_)                         ! J/m^3
    No2Si_V(UnitPoynting_)   = No2Si_V(UnitEnergyDens_)*No2Si_V(UnitU_)! J/m^2/s
    No2Si_V(UnitJ_)          = No2Si_V(UnitB_)/( No2Si_V(UnitX_)*cMu ) ! A/m^2
    No2Si_V(UnitElectric_)   = No2Si_V(UnitU_)*No2Si_V(UnitB_)         ! V/m
    No2Si_V(UnitTemperature_)= No2Si_V(UnitP_) &
         /( No2Si_V(UnitN_)*cBoltzmann )                               ! K 
    No2Si_V(UnitDivB_)       = No2Si_V(UnitB_)/No2Si_V(UnitX_)         ! T/m
    No2Si_V(UnitAngle_)      = 1.0                                     ! radian

    !\
    ! Set inverse conversion SI -> normalized
    !/
    Si2No_V = 1.0/No2Si_V

    ! As a default use SI units, so below only the differences need to be set
    Io2Si_V = 1.0
    No2Io_V = No2Si_V

    Io2Si_V(UnitX_)           = rSun                      ! R
    Io2Si_V(UnitRho_)         = 1.0E+3                    ! g/cm^3
    Io2Si_V(UnitN_)           = 1.0E+6                    ! #/cm^3
    Io2Si_V(UnitU_)           = 1.0E+3                    ! km/s
    Io2Si_V(UnitP_)           = 1.0E-1                    ! dyne/cm^2
    Io2Si_V(UnitB_)           = 1.0E-4                    ! Gauss
    Io2Si_V(UnitRhoU_)        = 1.0E+1                    ! g/cm^2/s
    Io2Si_V(UnitEnergydens_)  = 1.0E-1                    ! erg/cm^3
    Io2Si_V(UnitJ_)           = 1.0E-6                    ! uA/m^2
    Io2Si_V(UnitDivB_)        = 1.0E-2                    ! Gauss/cm
    Io2Si_V(UnitAngle_)       = cRadToDeg                 ! degrees

    ! Calculate the remaining unit conversions
    Si2Io_V = 1/Io2Si_V
    No2Io_V = No2Si_V*Si2Io_V
    Io2No_V = 1/No2Io_V

    Gbody  = -cGravitation*mSun*(Si2No_V(UnitU_)**2 * Si2No_V(UnitX_))

  end subroutine EEE_initialize

  !============================================================================

  subroutine EEE_set_parameters(NameCommand)

    use ModReadParam,     ONLY: read_var
    use EEE_ModGL98,      ONLY: set_parameters_GL98
    use EEE_ModTD99,      ONLY: set_parameters_TD99
    use EEE_ModArch,      ONLY: set_parameters_arch
    use EEE_ModShearFlow, ONLY: set_parameters_shearflow
    use EEE_ModCms,       ONLY: set_parameters_cms
    use EEE_ModCommonVariables, ONLY: &
         UseCme, DoAddFluxRope, UseTD, UseGL, UseShearFLow, UseArch, &
         DoAddFluxRope, LongitudeCme, LatitudeCme, OrientationCme, UseCms

    character(len=*), intent(in) :: NameCommand

    character(len=20):: TypeCme

    character(len=*), parameter:: NameSub = 'EEE_ModMain::EEE_set_parameters'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case("#CME")
       call read_var("UseCme", UseCme)
       if(UseCme)then
          call read_var("DoAddFluxRope",  DoAddFluxRope)
          call read_var("LongitudeCme",   LongitudeCme)
          call read_var("LatitudeCme",    LatitudeCme)
          call read_var("OrientationCme", OrientationCme)
          call read_var("TypeCme", TypeCme, IsUpperCase=.true.)
          select case(TypeCme)
          case("TITOV-DEMOULIN", "TD")
             UseTD = .true.
             call set_parameters_TD99(NameCommand)
          case("GIBSON-LOW", "GL")
             UseGL = .true.
             call set_parameters_GL98(NameCommand)
          case("BREAKOUT")
             UseShearFlow = .true.
             call set_parameters_shearflow(NameCommand)
             UseArch = .true.
             call set_parameters_arch(NameCommand)
          case default
             call CON_stop(NameSub//': invalid value for TypeCme='//TypeCme)
          end select
       end if

       ! The remaining commands are preserved for backwards compatibility
       ! and possibly for expert use (more options than #CME command)
    case("#ARCH")
       call set_parameters_arch(NameCommand)
    case("#TD99FLUXROPE")
       UseCme = .true.
       UseTD  = .true.
       call set_parameters_TD99(NameCommand)
    case("#GL98FLUXROPE")
       UseCme = .true.
       UseGL  = .true.
       call set_parameters_GL98(NameCommand)
    case("#SHEARFLOW")
       UseCme       = .true.
       UseShearFlow = .true.
       call set_parameters_shearflow(NameCommand)
    case("#CMS")
       UseCme = .true.
       UseCms = .true.     
       call set_parameters_cms(NameCommand)
    end select

  end subroutine EEE_set_parameters

  !============================================================================

  subroutine EEE_get_state_BC(x_D,Rho,U_D,B_D,p,Time,n_step,iteration_number)

    use EEE_ModCommonVariables, ONLY: UseCme, UseTD, UseShearFlow, UseGL, &
         UseCms
    use EEE_ModTD99, ONLY: get_transformed_TD99fluxrope, DoBqField
    use EEE_ModShearFlow, ONLY: get_shearflow
    use EEE_ModGL98, ONLY: get_GL98_fluxrope, adjust_GL98_fluxrope
    use EEE_ModCms, ONLY: get_cms

    real, intent(in) :: x_D(3), Time
    real, intent(out) :: Rho, U_D(3), B_D(3), p
    integer, intent(in) :: n_step,iteration_number

    real :: Rho1, U1_D(3), B1_D(3), p1
    !--------------------------------------------------------------------------

    ! initialize perturbed state variables
    Rho = 0.0; U_D = 0.0; B_D = 0.0; p = 0.0

    if(.not.UseCme) RETURN

    if (UseTD) then
       call get_transformed_TD99fluxrope(x_D,B1_D,&
            U1_D,n_step,Iteration_Number,Rho1,Time)

       if(.not.DoBqField) U1_D=0.0

       Rho = Rho + Rho1; U_D = U_D + U1_D; B_D = B_D + B1_D
    end if

    if(UseGL)then
       ! Add Gibson & Low (GL98) flux rope
       call get_GL98_fluxrope(x_D, Rho1, p1, B1_D)
       B_D = B_D + B1_D
    endif

    if(UseShearFlow)then
       call get_shearflow(x_D, Time, U1_D, iteration_number)

       U_D = U_D + U1_D
    end if

    if(UseCms) call get_cms(x_D, B_D)

  end subroutine EEE_get_state_BC

  !============================================================================

  subroutine EEE_get_state_init(x_D, Rho, B_D, p, n_step, iteration_number)

    use EEE_ModCommonVariables, ONLY: UseCme, DoAddFluxRope, UseTD, UseGL, &
         UseCms
    use EEE_ModGL98, ONLY: get_GL98_fluxrope, adjust_GL98_fluxrope
    use EEE_ModTD99, ONLY: get_transformed_TD99fluxrope
    use EEE_ModCms,  ONLY: get_cms

    real, intent(in) :: x_D(3)
    real, intent(out) :: Rho, B_D(3), p
    integer, intent(in) :: n_step,iteration_number

    real :: U_D(3)
    real :: Rho1, U1_D(3), B1_D(3), p1
    !--------------------------------------------------------------------------

    ! initialize perturbed state variables
    Rho = 0.0; U_D = 0.0; B_D = 0.0; p = 0.0

    if(.not. (UseCme .and. (DoAddFluxRope .or. UseCms))) RETURN

    if(UseTD)then
       ! Add Titov & Demoulin (TD99) flux rope
       call get_transformed_TD99fluxrope(x_D, B1_D, &
            U1_D, n_step, iteration_number, Rho1)
       Rho = Rho + Rho1; U_D = U_D + U1_D; B_D = B_D + B1_D
    endif

    if(UseGL)then
       ! Add Gibson & Low (GL98) flux rope
       call get_GL98_fluxrope(x_D, Rho1, p1, B1_D)
       call adjust_GL98_fluxrope(Rho1, p1)
       Rho = Rho + Rho1; B_D = B_D + B1_D; p = p + p1
    end if

    if(UseCms) call get_cms(x_D, B_D)

  end subroutine EEE_get_state_init

end Module EEE_ModMain
