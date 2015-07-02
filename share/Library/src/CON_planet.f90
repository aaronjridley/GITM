!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!
!BOP
!
!MODULE: CON_planet - planet parameters shared by components
!INTERFACE:
module CON_planet

  !DESCRIPTION:
  ! Physical information about the planet. The planet is described 
  ! with its name. Default values can be set with {\bf planet\_init}.
  ! Simplifying assumptions, such as no rotation, aligned magnetic
  ! and rotational axes etc. can be made.
  !
  ! This is a public class. The variables should be modified by CON only.
  ! Components can only access the data through the inquiry methods 
  ! via the {\bf CON\_physics} class.

  !USES:
  use ModNumConst, ONLY: cTwoPi
  use ModPlanetConst
  use ModTimeConvert, ONLY: TimeType, time_int_to_real

  !REVISION HISTORY:
  ! 01Aug03 - Aaron Ridly <ridley@umich.edu> and 
  !           Gabor Toth <gtoth@umich.edu>   - initial prototype/prolog/code
  ! 23Mar03 - added get_planet subroutine for OO type access
  ! 06May04 - K.C. Hansen and G. Toth added Saturn
  !           G.Toth fixed bugs in degree to radian conversions
  !EOP

  implicit none

  save

  character(len=*), parameter, private :: NameMod='CON_planet'

  character (len=lNamePlanet) :: NamePlanet = ''

  ! Define variables
  real           :: RadiusPlanet
  real           :: MassPlanet
  real           :: TiltRotation
  real           :: IonosphereHeight
  real           :: OmegaPlanet     ! Rotation + Orbit
  real           :: OmegaRotation    ! Rotation
  real           :: RotPeriodPlanet ! Rotation
  real           :: OmegaOrbit      ! Orbit
  real           :: AngleEquinox
  type(TimeType) :: TimeEquinox = TimeType(2000, 1, 1, 0, 0, 0, &
       0.0, 0.0_REAL8_, '20000101000000')

  ! Magnetic field type and strength in teslas
  character (len=lTypeBField) :: TypeBField = 'DIPOLE'
  real                        :: DipoleStrength
  real    :: MagAxisThetaGeo  ! Permanent theta  in GEO
  real    :: MagAxisPhiGeo    ! Permanent phi    in GEO

  ! Orientation of the axes
  real    :: RotAxisTheta      ! Permanent theta angle in GSE
  real    :: RotAxisPhi        ! Permanent phi   angle in GSE
  real    :: MagAxisTheta      ! Current   theta angla in GSE
  real    :: MagAxisPhi        ! Current   phi   angla in GSE

  ! Offset of the magnetic field center
  real    :: MagCenter_D(3) = (/ 0.0, 0.0, 0.0/)

  ! Optional changes relative to the "real" planet
  logical :: UseRotation     = .true.
  logical :: UseAlignedAxes  = .false.
  logical :: UseRealRotAxis  = .true.
  logical :: UseSetRotAxis   = .false.
  logical :: UseRealMagAxis  = .true.
  logical :: UseSetMagAxis   = .false.

  ! Frequency of updating the magnetic field information
  logical :: DoUpdateB0      = .true.
  real    :: DtUpdateB0      = 0.0001

  ! A primary axis is set to the true value
  ! a secondary axis is aligned with the primary axis
  logical :: IsRotAxisPrimary = .true., IsMagAxisPrimary = .true.

  ! A logical to indicate it the default parameters have been
  ! modified for this planet
  logical :: IsPlanetModified = .false.

contains

  !BOP ========================================================================
  !IROUTINE: set_planet_defaults - set Earth to be defaults
  !INTERFACE:
  subroutine set_planet_defaults

    !DESCRIPTION:
    ! Initialize parameters for Earth as the default planet.  This is in case
    ! there is no \#PLANET command.
    !EOP
    character (len=*), parameter :: NameSub = NameMod//'::set_planet_defaults'

    !-------------------------------------------------------------------------
    Planet_          = Earth_
    NamePlanet       = NamePlanet_I(Earth_)
    RadiusPlanet     = rPlanet_I(Earth_)
    MassPlanet       = mPlanet_I(Earth_)
    RotPeriodPlanet  = RotationPeriodPlanet_I(Earth_)
    TiltRotation     = TiltPlanet_I(Earth_)
    IonosphereHeight = IonoHeightPlanet_I(Earth_)
    OmegaOrbit       = cTwoPi/OrbitalPeriodPlanet_I(Earth_)
    OmegaPlanet      = OmegaOrbit + cTwoPi/RotationPeriodPlanet_I(Earth_)
    AngleEquinox     = cTwoPi * &
         ( iHourEquinoxPlanet_I(Earth_) * 3600 &
         + iMinuteEquinoxPlanet_I(Earth_) * 60 &
         + iSecondEquinoxPlanet_I(Earth_) &
         + FracSecondEquinoxPlanet_I(Earth_) &
         ) / (24 * 3600)
    TimeEquinox      = TimeType(&
         iYearEquinoxPlanet_I(Earth_), &
         iMonthEquinoxPlanet_I(Earth_), &
         iDayEquinoxPlanet_I(Earth_), &
         iHourEquinoxPlanet_I(Earth_),   &
         iMinuteEquinoxPlanet_I(Earth_), &
         iSecondEquinoxPlanet_I(Earth_), &
         FracSecondEquinoxPlanet_I(Earth_), &
         0.0_Real8_, '20000320073500')
    call time_int_to_real(TimeEquinox)
    TypeBField       = TypeBFieldPlanet_I(Earth_)
    DipoleStrength   = DipoleStrengthPlanet_I(Earth_)
    MagAxisThetaGeo  = bAxisThetaPlanet_I(Earth_)     
    MagAxisPhiGeo    = bAxisPhiPlanet_I(Earth_)       

    
  end subroutine set_planet_defaults

  !BOP ========================================================================
  !IROUTINE: is_planet_init - initialize parameters if planet is known
  !INTERFACE:
  function is_planet_init(NamePlanetIn) result(IsKnown)

    !INPUT ARGUMENTS:
    character(len=*), intent(in) :: NamePlanetIn
    
    !RETURN VALUE:
    logical :: IsKnown

    !DESCRIPTION:
    ! Initialize parameters for the planet identified by its name and
    ! return true if the planet is known. If it is not known return false. 
    ! Store the name in either case. The planet data can be initialized at most
    ! once.
    !EOP
    character (len=*), parameter :: NameSub = NameMod//'::is_planet_init'
    integer :: i

    logical :: IsInitialized = .false.
    !-------------------------------------------------------------------------

    IsKnown = .true.
    if(IsInitialized)then
       if(NamePlanet == NamePlanetIn) RETURN
       call CON_stop(NameSub//&
            ' ERROR: attempt to change planet name from '// &
            trim(NamePlanet)//' to '//NamePlanetIn)
    end if

    NamePlanet    = NamePlanetIn
    IsInitialized = .true.

    IsKnown       = .false.
    do i = NoPlanet_, MaxPlanet
      if (NamePlanet == NamePlanet_I(i)) then
         IsKnown = .true.
         Planet_ = i
         EXIT
      end if
    end do

    if (.not. IsKnown)  then
       Planet_    = NewPlanet_
       NamePlanet = NamePlanetIN
    end if

    ! Set all values for the selected planet
    RadiusPlanet     = rPlanet_I(Planet_)
    MassPlanet       = mPlanet_I(Planet_)
    RotPeriodPlanet  = RotationPeriodPlanet_I(Planet_)
    TiltRotation     = TiltPlanet_I(Planet_)
    IonosphereHeight = IonoHeightPlanet_I(Planet_)
    if (RotationPeriodPlanet_I(Planet_) == 0.0) then
       OmegaRotation = 0.0 
    else
       OmegaRotation = cTwoPi/RotationPeriodPlanet_I(Planet_)
    end if
    if (OrbitalPeriodPlanet_I(Planet_) == 0.0) then
       OmegaOrbit = 0.0 
    else
       OmegaOrbit = cTwoPi/OrbitalPeriodPlanet_I(Planet_)
    end if
    OmegaPlanet     = OmegaRotation + OmegaOrbit
    AngleEquinox =  &
         cTwoPi * ( iHourEquinoxPlanet_I(Planet_) * 3600 + iMinuteEquinoxPlanet_I(Planet_) * 60 &
         + iSecondEquinoxPlanet_I(Planet_) + FracSecondEquinoxPlanet_I(Planet_)) / (24 * 3600)
    TimeEquinox  = TimeType(&
         iYearEquinoxPlanet_I(Planet_), &
         iMonthEquinoxPlanet_I(Planet_), &
         iDayEquinoxPlanet_I(Planet_), &
         iHourEquinoxPlanet_I(Planet_), &
         iMinuteEquinoxPlanet_I(Planet_), &
         iSecondEquinoxPlanet_I(Planet_), &
         FracSecondEquinoxPlanet_I(Planet_), &
         0.0_Real8_, '')

    ! Magnetic field type and strength in teslas
    TypeBField        = TypeBFieldPlanet_I(Planet_)
    DipoleStrength    = DipoleStrengthPlanet_I(Planet_)
    MagAxisThetaGeo   = bAxisThetaPlanet_I(Planet_)  ! Permanent theta  in GEO
    MagAxisPhiGeo     = bAxisPhiPlanet_I(Planet_)    ! Permanent phi    in GEO

    ! For Enceladus the dipole is at Saturn's center
    if(Planet_==Enceladus_) MagCenter_D(2)    = 944.23

  end function is_planet_init

  !===========================================================================

  subroutine read_planet_var(NameCommand)

    use ModUtilities, ONLY: upper_case
    use ModReadParam, ONLY: read_var, lStringLine

    character (len=*), intent(in) :: NameCommand

    character (len=*), parameter :: NameSub = NameMod//'::read_planet_var'

    ! Planet related temporary variables
    character (len=lNamePlanet) :: NamePlanetIn
    character (len=lStringLine) :: NamePlanetCommands=''
    logical :: UseNonDipole
    !-------------------------------------------------------------------------

    select case(NameCommand)
    case("#PLANET", "#MOON", "#COMET")

       call read_var('NamePlanet',NamePlanetIn)
       call upper_case(NamePlanetIn)

       ! Check if planet has been already initialized
       if(NamePlanet == NamePlanetIn) RETURN

       if (NamePlanetCommands /= '') &
            call CON_stop(NameSub// &
            ' ERROR: #PLANET should precede '// &
            NamePlanetCommands)

       if ( .not. is_planet_init(NamePlanetIn) ) then

          call read_var('RadiusPlanet', RadiusPlanet)
          call read_var('MassPlanet',   MassPlanet)
          call read_var('OmegaPlanet',  OmegaPlanet)
          if (OmegaPlanet /= 0.0) then
             RotPeriodPlanet = cTwoPi/OmegaPlanet
          else
             RotPeriodPlanet = 0.0
          end if
          call read_var('TiltRotation', TiltRotation)
          TiltRotation = TiltRotation * cDegToRad
          call read_var('TypeBField',   TypeBField)

          call upper_case(TypeBField)

          select case(TypeBField)

          case('NONE')

             MagAxisTheta   = 0.0
             MagAxisPhi     = 0.0
             UseSetMagAxis  = .true.
             UseRealMagAxis = .false.
             DipoleStrength = 0.0

          case('DIPOLE','QUADRUPOLE','OCTUPOLE')

             call read_var('MagAxisThetaGeo', MagAxisThetaGeo)
             MagAxisThetaGeo = MagAxisThetaGeo * cDegToRad
             call read_var('MagAxisPhiGeo',   MagAxisPhiGeo)
             MagAxisPhiGeo = MagAxisPhiGeo * cDegToRad
             call read_var('DipoleStrength',DipoleStrength)

             if (TypeBField == 'QUADRUPOLE') then
                call CON_stop(NameSub// &
                     ' ERROR: quadrupole field unimplemented')
             endif

             if (TypeBField == 'OCTUPOLE') then
                call CON_stop(NameSub// &
                     ' ERROR: octupole field unimplemented')
             endif

          case default
             call CON_stop(NameSub// &
                  ' ERROR: unimplemented TypeBField='//TypeBField)

          end select

       endif

    case('#IDEALAXES')
       ! This is a short version of setting one axis parallel with Z
       ! and the other one aligned with it

       NamePlanetCommands = '#IDEALAXES ' // NamePlanetCommands
       IsPlanetModified = .true.

       UseRealRotAxis   = .false.
       IsRotAxisPrimary = .true.
       UseSetRotAxis    = .true.
       RotAxisTheta     = 0.0
       RotAxisPhi       = 0.0
       UseRealMagAxis   = .false.
       IsMagAxisPrimary = .false.
       UseSetMagAxis    = .false.

    case('#ROTATIONAXIS')

       NamePlanetCommands = '#ROTATIONAXIS ' // NamePlanetCommands
       IsPlanetModified = .true.
       UseRealRotAxis = .false.

       call read_var('IsRotAxisPrimary', IsRotAxisPrimary)
       if (IsRotAxisPrimary) then

          UseSetRotAxis = .true.

          call read_var('RotAxisTheta', RotAxisTheta)
          if(RotAxisTheta < 0.0)call CON_stop(NameSub// &
               ' ERROR: negative tilt should be entered as Phi=180.0')
          RotAxisTheta = cDegToRad * RotAxisTheta

          call read_var('RotAxisPhi', RotAxisPhi)
          RotAxisPhi = cDegToRad * RotAxisPhi
       else
          if(.not.IsMagAxisPrimary)call CON_stop(NameSub// &
               ' ERROR: either rotation or magnetic axis must be primary')
       end if

    case('#MAGNETICAXIS')

       NamePlanetCommands = '#MAGNETICAXIS ' // NamePlanetCommands
       IsPlanetModified = .true.
       UseRealMagAxis = .false.

       call read_var('IsMagAxisPrimary', IsMagAxisPrimary)
       if (IsMagAxisPrimary) then

          UseSetMagAxis = .true.

          call read_var('MagAxisTheta', MagAxisTheta)
          if(MagAxisTheta < 0.0)call CON_stop(NameSub// &
               ' ERROR: negative tilt should be entered as Phi=180.0')
          MagAxisTheta = cDegToRad * MagAxisTheta

          call read_var('MagAxisPhi', MagAxisPhi)
          MagAxisPhi = cDegToRad * MagAxisPhi
       else
          if(.not.IsRotAxisPrimary)call CON_stop(NameSub//&
               ' ERROR: either rotation or magnetic axis must be primary')
       end if

    case('#MAGNETICCENTER')

       NamePlanetCommands = '#MAGNETICCENTER ' // NamePlanetCommands
       IsPlanetModified = .true.
       call read_var('MagneticCenterX',MagCenter_D(1))
       call read_var('MagneticCenterY',MagCenter_D(2))
       call read_var('MagneticCenterZ',MagCenter_D(3))

    case('#ROTATION')

       NamePlanetCommands = '#ROTATION ' // NamePlanetCommands
       IsPlanetModified = .true.

       call read_var('UseRotation', UseRotation)
       if (.not.UseRotation) then
          OmegaPlanet     = 0.0
          RotPeriodPlanet = 0.0
       else
          call read_var('Rotation period [hours]',  RotPeriodPlanet)
          RotPeriodPlanet = RotPeriodPlanet * 3600
          OmegaPlanet = cTwoPi / RotPeriodPlanet
       endif

    case('#NONDIPOLE')

       NamePlanetCommands = '#NONDIPOLE ' // NamePlanetCommands
       IsPlanetModified = .true.

       call read_var('UseNonDipole',UseNonDipole)
       if (.not.UseNonDipole) then
          TypeBField = 'DIPOLE'
       else
          call CON_stop(NameSub// &
               ' ERROR: nondipole magnetic field unimplemented')
       endif

    case('#DIPOLE')

       NamePlanetCommands = '#DIPOLE ' // NamePlanetCommands
       IsPlanetModified = .true.

       call read_var('DipoleStrength',DipoleStrength)
       if(DipoleStrength ==0.0)then
          TypeBField="NONE"
       else
          TypeBField="DIPOLE"
       end if
       
    case('#UPDATEB0')

       call read_var('DtUpdateB0',DtUpdateB0)

    end select

  end subroutine read_planet_var

  !==========================================================================

  subroutine check_planet_var(IsProc0, DoTimeAccurate)

    logical, intent(in) :: IsProc0, DoTimeAccurate

    character (len=*), parameter :: NameSub=NameMod//'::check_planet_var'

    ! The rotation and magnetic axes are aligned if any of them is not a 
    ! primary axis.
    UseAlignedAxes = (.not. IsRotAxisPrimary) .or. (.not. IsMagAxisPrimary)

    ! Warn if setting is unphysical
    if(UseSetMagAxis .and. UseRealRotAxis .and. IsProc0) &
         write(*,*)NameSub,' WARNING: magnetic axis is explicitly set ',&
         'while rotation axis is calculated from real time ?!'

    ! Check if there is a need to update the magnetic field
    DoUpdateB0 = DtUpdateB0 > 0.0 .and. DoTimeAccurate .and. UseRotation &
         .and. .not.UseAlignedAxes

  end subroutine check_planet_var

  !==========================================================================

  subroutine get_planet( &
       NamePlanetOut, RadiusPlanetOut, MassPlanetOut, OmegaPlanetOut, &
       RotationPeriodOut, IonosphereHeightOut, &
       UseRotationOut, DipoleStrengthOut, DoUpdateB0Out, DtUpdateB0Out)

    character(len=*), optional, intent(out) :: NamePlanetOut
    real,             optional, intent(out) :: RadiusPlanetOut
    real,             optional, intent(out) :: MassPlanetOut
    real,             optional, intent(out) :: OmegaPlanetOut
    real,             optional, intent(out) :: RotationPeriodOut
    real,             optional, intent(out) :: IonosphereHeightOut
    logical,          optional, intent(out) :: UseRotationOut
    real,             optional, intent(out) :: DipoleStrengthOut
    logical,          optional, intent(out) :: DoUpdateB0Out
    real,             optional, intent(out) :: DtUpdateB0Out
    !-----------------------------------------------------------------------
    if(present(NamePlanetOut))      NamePlanetOut       = NamePlanet
    if(present(RadiusPlanetOut))    RadiusPlanetOut     = RadiusPlanet
    if(present(MassPlanetOut))      MassPlanetOut       = MassPlanet
    if(present(OmegaPlanetOut))     OmegaPlanetOut      = OmegaPlanet
    if(present(RotationPeriodOut))  RotationPeriodOut   = RotPeriodPlanet
    if(present(IonosphereHeightOut))IonosphereHeightOut = IonosphereHeight
    if(present(UseRotationOut))     UseRotationOut      = UseRotation
    if(present(DipoleStrengthOut))  DipoleStrengthOut   = DipoleStrength
    if(present(DoUpdateB0Out))      DoUpdateB0Out       = DoUpdateB0
    if(present(DtUpdateB0Out))      DtUpdateB0Out       = DtUpdateB0

  end subroutine get_planet

end module CON_planet
